#!/usr/bin/env python3
"""Independently recompute and validate the TF48 Stage 4 package."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


COMPARTMENTS = ["tumor_like", "immune", "stroma_fibroblast", "endothelial", "ambiguous_or_low_signal"]
KEY_FIELDS = ["dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment", "ligand", "receptor", "pathway"]
RECURRENCE_KEY = ["dataset", "sender_compartment", "receiver_compartment", "ligand", "receptor"]
F5_ROW_KEY = ["sample_id", "replicate", "condition_id", "axis_id", "sender_compartment", "receiver_compartment", "ligand", "receptor"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def close(left: object, right: object, tolerance: float = 2e-6) -> bool:
    return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=2e-8)


def same_construct_value(field: str, left: object, right: object) -> bool:
    try:
        if field in {"sender_cells", "receiver_cells"}:
            return int(left) == int(right)
        if field in {"ligand_addition", "receptor_addition", "reporter_addition"}:
            left_value = float(left)
            right_value = float(right)
            return (
                math.isfinite(left_value)
                and math.isfinite(right_value)
                and math.isclose(left_value, right_value, rel_tol=1e-14, abs_tol=1e-14)
            )
    except (TypeError, ValueError):
        return False
    return str(left) == str(right)


def pilot_tier(row: dict[str, str], score: float) -> str:
    null_ok = all(row.get(field, "") and float(row[field]) <= 0.10 for field in ["spatial_null_p", "label_null_p", "fake_lr_null_p"])
    if score >= 0.75 and null_ok and int(row["external_method_support_count"]) == 2:
        return "high_pilot_v3"
    if score >= 0.55:
        return "medium_pilot_v3"
    return "low_pilot_v3"


def target_metrics_from_stats(
    genes: list[str], receiver: str, sample_stats: dict[tuple[str, str], dict[str, float]],
) -> dict[str, object]:
    detected = [gene for gene in genes if sample_stats.get((gene, COMPARTMENTS[0]), {}).get("gene_detected_in_matrix", 0)]
    means = {}
    positive_fractions = {}
    for comp in COMPARTMENTS:
        spots = int(sample_stats[(genes[0], comp)]["compartment_spots"])
        total = sum(float(sample_stats[(gene, comp)]["gene_expression_sum"]) for gene in detected)
        means[comp] = total / max(spots * max(len(detected), 1), 1)
        positive_union = 0
        for gene in detected:
            positive_union |= int(sample_stats[(gene, comp)]["positive_spot_bitset_hex"] or "0", 16)
        positive_fractions[comp] = positive_union.bit_count() / spots if spots else 0.0
    other_spots = sum(int(sample_stats[(genes[0], comp)]["compartment_spots"]) for comp in COMPARTMENTS if comp != receiver)
    other_total = sum(float(sample_stats[(gene, comp)]["gene_expression_sum"]) for gene in detected for comp in COMPARTMENTS if comp != receiver)
    other_mean = other_total / max(other_spots * max(len(detected), 1), 1)
    receiver_mean = means[receiver]
    percentile = sum(value <= receiver_mean for value in means.values()) / len(COMPARTMENTS)
    enrichment = (receiver_mean + 1e-9) / (other_mean + 1e-9)
    support_without_positive_fraction = int(len(detected) >= 2 and percentile >= 0.75 and enrichment >= 1.25)
    return {
        "target_genes_in_panel": len(genes),
        "target_genes_detected": len(detected),
        "target_genes_detected_list": ",".join(detected),
        "receiver_target_mean": receiver_mean,
        "other_compartment_target_mean": other_mean,
        "receiver_target_enrichment": enrichment,
        "receiver_target_percentile": percentile,
        "target_support": support_without_positive_fraction,
        "receiver_target_positive_spot_fraction": positive_fractions[receiver],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--source-lock", required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    project = Path(args.project).resolve()
    root = Path(args.result_dir).resolve()
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({"check": name, "passed": int(bool(passed)), "observed": observed, "expected": expected})

    summary = json.loads((root / "summary.json").read_text())
    check("builder_status", summary["status"] == "PASS_TF48_STAGE4_REVISION_ANALYSES_BUILT", summary["status"], "PASS_TF48_STAGE4_REVISION_ANALYSES_BUILT")
    check("source_commit", summary["source_commit"] == args.source_commit, summary["source_commit"], args.source_commit)
    check("failures_empty", len(read_tsv(root / "failures.tsv")) == 0, len(read_tsv(root / "failures.tsv")), 0)

    lock = read_tsv(Path(args.source_lock))
    lock_errors = []
    for row in lock:
        if row["role"] == "git_head":
            observed = subprocess.check_output(["git", "-C", str(project), "rev-parse", "HEAD"], text=True).strip()
        else:
            path = project / row["path"]
            observed = digest(path) if path.is_file() else "MISSING"
        if observed != row["sha256"]:
            lock_errors.append(row["path"])
    check("source_lock_nonempty", len(lock) >= 49, len(lock), ">=49")
    check("source_lock_rehash", not lock_errors, ";".join(lock_errors), "all rows match live committed inputs")

    construct = read_tsv(root / "semisim_construct_provenance.tsv")
    construct_keys = [tuple(row[field] for field in F5_ROW_KEY) for row in construct]
    check("construct_rows_unique", len(construct) == 1140 and len(construct_keys) == len(set(construct_keys)), f"{len(construct)}/{len(set(construct_keys))}", "1140/1140")
    manifest_hashes = {}
    truth_rows = {}
    batches = project / "results/task_f_cross_platform/f5_semisim_full_input_batches"
    for batch in sorted(path for path in batches.iterdir() if path.is_dir()):
        for row in read_tsv(batch / "generated_input_manifest.tsv"):
            manifest_hashes[(row["sample_id"], row["replicate"], row["condition_id"], row["role"])] = row["sha256"]
        for row in read_tsv(batch / "truth_contract.tsv"):
            truth_rows[tuple(row[field] for field in F5_ROW_KEY)] = row
    construct_mismatch = 0
    actual_hash_mismatch = 0
    for row in construct:
        key = tuple(row[field] for field in F5_ROW_KEY)
        truth = truth_rows.get(key)
        if not truth:
            construct_mismatch += 1
            continue
        for field in ["truth_class", "resource_axis", "sender_cells", "receiver_cells", "ligand_addition", "receptor_addition", "reporter_addition"]:
            if not same_construct_value(field, row[field], truth[field]):
                construct_mismatch += 1
                break
        for role, field in [("shared_input.h5ad", "h5ad_sha256"), ("matrix.mtx.gz", "matrix_sha256"), ("features.tsv.gz", "features_sha256"), ("barcodes.tsv.gz", "barcodes_sha256"), ("meta.tsv.gz", "metadata_sha256")]:
            if row[field] != manifest_hashes.get((row["sample_id"], row["replicate"], row["condition_id"], role)):
                construct_mismatch += 1
                break
            referenced_path = project / row[field.replace("sha256", "path")]
            if not referenced_path.is_file() or digest(referenced_path) != row[field]:
                actual_hash_mismatch += 1
                break
    check("construct_recomputed_from_retained_manifests", construct_mismatch == 0 and len(truth_rows) == 1140, f"mismatch={construct_mismatch};truth={len(truth_rows)}", "mismatch=0;truth=1140")
    check("construct_actual_artifact_rehash", actual_hash_mismatch == 0, actual_hash_mismatch, 0)

    resolution = read_tsv(root / "finite_null_resolution.tsv")
    check("null_resolution", len(resolution) == 6 and sorted(int(r["draws"]) for r in resolution) == [30, 50, 100, 999, 999, 999], f"rows={len(resolution)};draws={sorted(int(r['draws']) for r in resolution)}", "6 rows; 30,50,100,999,999,999")
    nominal = next(row for row in read_tsv(root / "pilot_null_threshold_sensitivity.tsv") if row["nominal_threshold"] == "0.1")
    check("nominal_missing_fail_closed", nominal["missing_fake_rows_fail_closed"] == "9", nominal["missing_fake_rows_fail_closed"], 9)

    pilot_source = read_tsv(project / "scores/spatiallr_trust_score_pilot_v3.tsv")
    pilot_output = {tuple(row[field] for field in KEY_FIELDS): row for row in read_tsv(root / "pilot_v3_loo_recurrence.tsv")}
    sample_sets: dict[str, set[str]] = defaultdict(set)
    passing: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in pilot_source:
        sample_sets[row["dataset"]].add(row["sample_id"])
        if float(row["spatial_null_p"]) <= 0.10 or float(row["label_null_p"]) <= 0.10:
            passing[tuple(row[field] for field in RECURRENCE_KEY)].add(row["sample_id"])
    pilot_mismatch = 0
    pilot_high = pilot_changed = 0
    for row in pilot_source:
        key = tuple(row[field] for field in KEY_FIELDS)
        rec_key = tuple(row[field] for field in RECURRENCE_KEY)
        contributors = passing[rec_key]
        train_pass = len(contributors) - int(row["sample_id"] in contributors)
        train_total = len(sample_sets[row["dataset"]]) - 1
        support = min(1.0, 2 * train_pass / train_total)
        score = 0.40 * float(row["null_support_mean"]) + 0.15 * float(row["pilot_score_percentile"]) + 0.15 * support + 0.10 * float(row["annotation_quality"]) + 0.20 * float(row["external_method_support_fraction"])
        tier = pilot_tier(row, score)
        result = pilot_output.get(key)
        if not result or int(result["training_samples_passed"]) != train_pass or not close(result["loo_score"], score) or result["loo_tier"] != tier:
            pilot_mismatch += 1
        pilot_high += tier == "high_pilot_v3"
        pilot_changed += tier != row["confidence_tier_pilot_v3"]
    check("pilot_loo_independent_recompute", pilot_mismatch == 0 and len(pilot_output) == 4225, f"mismatch={pilot_mismatch};rows={len(pilot_output)}", "mismatch=0;rows=4225")
    check("pilot_loo_outcome", pilot_high == 217 and pilot_changed == 48, f"high={pilot_high};changed={pilot_changed}", "high=217;changed=48")
    check("pilot_units_no_technical_replicates", all(len(samples) in {4, 12, 18} for samples in sample_sets.values()), str({key: len(value) for key, value in sample_sets.items()}), "dataset sample counts 4,12,18")

    f5_source = read_tsv(project / "scores/spatiallr_trust_score_f5_semisim_full.tsv")
    common = [row for row in f5_source if row["resource_axis"] == "1" and row["truth_class"] in {"positive", "negative"} and row["method_scope_count"] == "3"]
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in common:
        groups[(row["sample_id"], row["replicate"])].append(row)
    check("f5_unit_structure", len(groups) == 15 and {len(value) for value in groups.values()} == {34}, f"groups={len(groups)};sizes={sorted(set(map(len, groups.values())))}", "15 groups; 34 rows each")
    scorers = {"COMMOT": "commot_score", "LIANA": "liana_score", "CellChat": "cellchat_score", "SpatialLR-Trust": "spatiallr_trust_score_f5_full"}
    expected_tp = {}
    for name, column in scorers.items():
        selected = []
        for group in groups.values():
            selected.extend(sorted(group, key=lambda row: (-float(row[column]), row["condition_id"], row["axis_id"]))[:4])
        expected_tp[name] = sum(row["truth_class"] == "positive" for row in selected)
    matched = {row["scorer"]: row for row in read_tsv(root / "f5_matched_operating_points.tsv") if row["sample_id"] == "all"}
    check("matched_operating_points_recomputed", expected_tp == {"COMMOT": 59, "LIANA": 56, "CellChat": 34, "SpatialLR-Trust": 46} and all(int(matched[name]["true_positive_rows"]) == count for name, count in expected_tp.items()), str(expected_tp), "COMMOT=59;LIANA=56;CellChat=34;SpatialLR-Trust=46")

    f5_output = {tuple(row[field] for field in F5_ROW_KEY): row for row in read_tsv(root / "f5_loo_recurrence.tsv")}
    f5_mismatch = f5_high = f5_changed = 0
    for row in f5_source:
        if row["resource_axis"] != "1":
            continue
        passed = {value for value in row["passed_biological_sample_ids"].split(";") if value}
        train_pass = int(row["joint_pair_null_pass_biological_samples"]) - int(row["sample_id"] in passed)
        fraction = train_pass / (int(row["biological_samples_total"]) - 1)
        score = 0.50 * float(row["pair_null_support_fraction"]) + 0.20 * float(row["method_support_fraction_all"]) + 0.20 * fraction + 0.10 * float(row["global_null_support_fraction"])
        high = score >= 0.75 and float(row["observed_local_lr_score"]) > 0 and int(row["pair_null_pass_count"]) == 3 and int(row["methods_detected"]) >= 2 and train_pass >= 2
        medium = not high and score >= 0.50 and float(row["observed_local_lr_score"]) > 0 and int(row["pair_null_pass_count"]) >= 2 and int(row["methods_detected"]) >= 1
        tier = "high_f5_full" if high else "medium_f5_full" if medium else "low_f5_full"
        result = f5_output.get(tuple(row[field] for field in F5_ROW_KEY))
        if not result or not close(result["loo_score"], score) or result["loo_tier"] != tier:
            f5_mismatch += 1
        f5_high += tier == "high_f5_full"
        f5_changed += tier != row["confidence_tier_f5_full"]
    check("f5_loo_independent_recompute", f5_mismatch == 0 and len(f5_output) == 1140, f"mismatch={f5_mismatch};rows={len(f5_output)}", "mismatch=0;rows=1140")
    check("f5_loo_outcome", f5_high == 41 and f5_changed == 40, f"high={f5_high};changed={f5_changed}", "high=41;changed=40")

    panel: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(project / "docs/pilot_pathway_target_panel.tsv"):
        gene = row["gene"].upper()
        if gene not in panel[row["pathway"]]:
            panel[row["pathway"]].append(gene)
    stats_rows = read_tsv(root / "target_gene_compartment_sufficient_statistics.tsv")
    expected_stats = 34 * len({gene for genes in panel.values() for gene in genes}) * len(COMPARTMENTS)
    check("target_sufficient_statistics_rows", len(stats_rows) == expected_stats, len(stats_rows), expected_stats)
    stats_by_sample: dict[tuple[str, str], dict[tuple[str, str], dict[str, float]]] = defaultdict(dict)
    for row in stats_rows:
        stats_by_sample[(row["dataset"], row["sample_id"])][(row["gene"], row["compartment"])] = {
            "compartment_spots": int(row["compartment_spots"]),
            "gene_expression_sum": float(row["gene_expression_sum"]),
            "gene_positive_spots": int(row["gene_positive_spots"]),
            "gene_detected_in_matrix": int(row["gene_detected_in_matrix"]),
            "positive_spot_bitset_hex": row["positive_spot_bitset_hex"],
        }
    target_output = {tuple(row[field] for field in KEY_FIELDS): row for row in read_tsv(root / "target_candidate_excluded.tsv")}
    current_target = {tuple(row[field] for field in KEY_FIELDS): row for row in read_tsv(project / "results/task_e_target_activation/target_activation_all_score_v3_with_target.tsv")}
    target_mismatch = original_source_mismatch = 0
    original_supported = high_supported = mif_supported = 0
    for source in pilot_source:
        key = tuple(source[field] for field in KEY_FIELDS)
        result = target_output.get(key)
        genes = panel[source["pathway"]]
        sample_stats = stats_by_sample[(source["dataset"], source["sample_id"])]
        original = target_metrics_from_stats(genes, source["receiver_compartment"], sample_stats)
        removed = {source["ligand"].upper(), source["receptor"].upper()}
        excluded = target_metrics_from_stats([gene for gene in genes if gene not in removed], source["receiver_compartment"], sample_stats)
        current = current_target[key]
        if int(current["target_activation_support"]) != original["target_support"] or not close(current["receiver_target_enrichment"], original["receiver_target_enrichment"]) or not close(current["receiver_target_positive_spot_fraction"], original["receiver_target_positive_spot_fraction"]):
            original_source_mismatch += 1
        if not result or int(result["excluded_target_support"]) != excluded["target_support"] or result["excluded_target_genes_detected_list"] != excluded["target_genes_detected_list"] or not close(result["excluded_receiver_target_enrichment"], excluded["receiver_target_enrichment"]) or not close(result["excluded_receiver_target_positive_spot_fraction"], excluded["receiver_target_positive_spot_fraction"]):
            target_mismatch += 1
        original_supported += original["target_support"]
        high_supported += original["target_support"] if source["confidence_tier_pilot_v3"] == "high_pilot_v3" else 0
        mif_supported += original["target_support"] if source["pathway"] == "MIF" else 0
    check("target_stats_reproduce_frozen_source", original_source_mismatch == 0 and original_supported == 1266 and high_supported == 134 and mif_supported == 107, f"mismatch={original_source_mismatch};all={original_supported};high={high_supported};mif={mif_supported}", "mismatch=0;all=1266;high=134;mif=107")
    check("target_exclusion_independent_recompute", target_mismatch == 0 and len(target_output) == 4225, f"mismatch={target_mismatch};rows={len(target_output)}", "mismatch=0;rows=4225")
    overlap = read_tsv(root / "target_panel_overlap.tsv")
    check("target_overlap", len(overlap) == 18 and any(row["pathway"] == "MIF" and row["overlap_gene"] == "CD74" for row in overlap), f"rows={len(overlap)};mif_cd74={any(row['pathway']=='MIF' and row['overlap_gene']=='CD74' for row in overlap)}", "rows=18;mif_cd74=True")
    mif_loo = read_tsv(root / "mif_target_gene_loo.tsv")
    check("featured_mif_gene_loo", sum(int(row["featured_case"]) for row in mif_loo) == 7, sum(int(row["featured_case"]) for row in mif_loo), 7)

    passed = sum(int(row["passed"]) for row in checks)
    failed = len(checks) - passed
    status = "PASS_TF48_STAGE4_REVISION_ANALYSES" if failed == 0 else "FAIL_TF48_STAGE4_REVISION_ANALYSES"
    with (root / "validation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "passed", "observed", "expected"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)
    validation = {"status": status, "checks": len(checks), "passed": passed, "failed": failed, "source_commit": args.source_commit, "authorized_next_action": "stage4_bilingual_manuscript_and_response_revision" if failed == 0 else "none"}
    (root / "validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    if failed:
        for row in checks:
            if not row["passed"]:
                print(json.dumps(row, sort_keys=True), file=sys.stderr)
        raise SystemExit(json.dumps(validation, sort_keys=True))
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
