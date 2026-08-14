#!/usr/bin/env python3
"""Build the prospective TF48 Stage 4 sensitivity and provenance package."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


COMPARTMENTS = [
    "tumor_like",
    "immune",
    "stroma_fibroblast",
    "endothelial",
    "ambiguous_or_low_signal",
]
KEY_FIELDS = [
    "dataset", "sample_id", "cancer", "sender_compartment",
    "receiver_compartment", "ligand", "receptor", "pathway",
]
RECURRENCE_KEY = [
    "dataset", "sender_compartment", "receiver_compartment", "ligand", "receptor",
]
F5_ROW_KEY = [
    "sample_id", "replicate", "condition_id", "axis_id",
    "sender_compartment", "receiver_compartment", "ligand", "receptor",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def as_float(value: object, default: float = 0.0) -> float:
    if value in {None, ""} or (isinstance(value, float) and math.isnan(value)):
        return default
    return float(value)


def roc_auc(labels: pd.Series, scores: pd.Series) -> float:
    labels_np = labels.astype(int).to_numpy()
    positive = labels_np == 1
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if not n_pos or not n_neg:
        return float("nan")
    ranks = scores.astype(float).rank(method="average").to_numpy()
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(labels: pd.Series, scores: pd.Series) -> float:
    frame = pd.DataFrame({"label": labels.astype(int), "score": scores.astype(float)})
    grouped = frame.groupby("score")["label"].agg(["sum", "count"]).sort_index(ascending=False)
    tp = grouped["sum"].cumsum()
    fp = (grouped["count"] - grouped["sum"]).cumsum()
    recall = tp / frame["label"].sum()
    precision = tp / (tp + fp)
    return float(((recall - recall.shift(fill_value=0)) * precision).sum())


def build_construct_provenance(project: Path, output: Path, source_commit: str) -> int:
    preflight = project / "results/task_f_cross_platform/f5_semisim_preflight"
    manifest = pd.read_csv(preflight / "full_manifest.tsv", sep="\t", dtype={"sample_id": str})
    inputs = pd.read_csv(preflight / "input_contract.tsv", sep="\t", dtype={"sample_id": str})
    pairs = pd.read_csv(preflight / "pair_truth_contract.tsv", sep="\t")
    input_hash = inputs.set_index("sample_id")["sha256"].to_dict()
    pair_fields = pairs.set_index("axis_id")[["synthetic_reporter", "fake_ligand", "fake_receptor", "selection_boundary"]].to_dict("index")
    generator = project / "scripts/server/generate_task_f_f5_semisim_batch.py"
    smoke_generator = project / "scripts/server/generate_task_f_f5_semisim_smoke.py"
    batches = project / "results/task_f_cross_platform/f5_semisim_full_input_batches"
    rows: list[dict[str, object]] = []
    for batch in sorted(path for path in batches.iterdir() if path.is_dir()):
        truth = pd.read_csv(batch / "truth_contract.tsv", sep="\t", dtype={"sample_id": str})
        files = pd.read_csv(batch / "generated_input_manifest.tsv", sep="\t", dtype={"sample_id": str})
        for file_row in files.itertuples(index=False):
            actual = project / file_row.path
            if not actual.is_file() or digest(actual) != file_row.sha256:
                raise ValueError(f"Generated F5 artifact hash mismatch: {file_row.path}")
        file_map = {
            (str(row.sample_id), int(row.replicate), row.condition_id, row.role): (row.path, row.sha256)
            for row in files.itertuples(index=False)
        }
        joined = truth.merge(
            manifest,
            on=["sample_id", "replicate", "condition_id", "condition_family", "truth_class"],
            how="left",
            validate="many_to_one",
        )
        for row in joined.to_dict("records"):
            identity = (str(row["sample_id"]), int(row["replicate"]), row["condition_id"])
            pair = pair_fields[row["axis_id"]]
            out = dict(row)
            out.update({
                "input_sha256": input_hash[str(row["sample_id"])],
                "generator_script": "scripts/server/generate_task_f_f5_semisim_batch.py",
                "generator_source_commit": source_commit,
                "generator_sha256": digest(generator),
                "region_definition": "75-cell nearest cross-compartment regions; x-axis-separated regions for the spatial hard negative",
                "signal_definition": "strength*max(base_gene_population_sd,1); ligand in sender; receptor and reporter in receiver",
                "operation_order": "signal injection; coordinate mixing; label mixing; expression-profile mixing; serialization",
                "synthetic_reporter": pair["synthetic_reporter"],
                "fake_ligand": pair["fake_ligand"],
                "fake_receptor": pair["fake_receptor"],
                "selection_boundary": pair["selection_boundary"],
                "smoke_generator_sha256": digest(smoke_generator),
            })
            for role, prefix in [
                ("shared_input.h5ad", "h5ad"), ("matrix.mtx.gz", "matrix"),
                ("features.tsv.gz", "features"), ("barcodes.tsv.gz", "barcodes"),
                ("meta.tsv.gz", "metadata"),
            ]:
                path_value, hash_value = file_map[(*identity, role)]
                out[f"{prefix}_path"] = path_value
                out[f"{prefix}_sha256"] = hash_value
            rows.append(out)
    if len(rows) != 1140:
        raise ValueError(f"Expected 1,140 construct rows, found {len(rows)}")
    fields = [
        "sample_id", "replicate", "condition_id", "condition_family", "truth_class", "axis_id",
        "sender_compartment", "receiver_compartment", "ligand", "receptor", "synthetic_reporter",
        "fake_ligand", "fake_receptor", "resource_axis", "seed", "spike_strength_sd",
        "coordinate_mix_fraction", "label_mix_fraction", "expression_mix_fraction",
        "sender_cells", "receiver_cells", "ligand_addition", "receptor_addition", "reporter_addition",
        "median_cross_distance_before_corruption", "input_sha256", "generator_script",
        "generator_source_commit", "generator_sha256", "smoke_generator_sha256", "region_definition",
        "signal_definition", "operation_order", "selection_boundary", "h5ad_path", "h5ad_sha256",
        "matrix_path", "matrix_sha256", "features_path", "features_sha256", "barcodes_path",
        "barcodes_sha256", "metadata_path", "metadata_sha256",
    ]
    write_tsv(output / "semisim_construct_provenance.tsv", rows, fields)
    return len(rows)


def pilot_tier(row: dict[str, str], score: float, threshold: float) -> str:
    null_ok = all(row.get(field, "") and float(row[field]) <= threshold for field in [
        "spatial_null_p", "label_null_p", "fake_lr_null_p",
    ])
    both = int(row["external_method_support_count"]) == 2
    if score >= 0.75 and null_ok and both:
        return "high_pilot_v3"
    if score >= 0.55:
        return "medium_pilot_v3"
    return "low_pilot_v3"


def build_null_and_pilot_loo(project: Path, output: Path) -> dict[str, int]:
    score_rows = read_tsv(project / "scores/spatiallr_trust_score_pilot_v3.tsv")
    score_keys = [tuple(row[field] for field in KEY_FIELDS) for row in score_rows]
    if len(score_keys) != len(set(score_keys)):
        raise ValueError("Pilot score keys are not unique")
    resolution = []
    specs = [
        ("pilot_v3", "spatial_coordinate", 50, 0, "none; candidate-level pilot screen"),
        ("pilot_v3", "compartment_label", 30, 0, "none; candidate-level pilot screen"),
        ("pilot_v3", "expression_matched_fake_lr", 100, 9, "none; 9 candidates lack a usable matched background"),
        ("f5_semisim", "spatial_coordinate", 999, 0, "BH within 49 directions per LR; secondary global BH across 833 rows"),
        ("f5_semisim", "compartment_label", 999, 0, "BH within 49 directions per LR; secondary global BH across 833 rows"),
        ("f5_semisim", "expression_matched_fake_lr", 999, 0, "BH within 49 directions per LR; secondary global BH across 833 rows"),
    ]
    for layer, model, draws, missing, correction in specs:
        denominator = draws + 1
        resolution.append({
            "analysis_layer": layer,
            "null_model": model,
            "draws": draws,
            "plus_one_formula": f"(exceedances+1)/({draws}+1)",
            "minimum_attainable_p": 1 / denominator,
            "grid_step": 1 / denominator,
            "missing_background_rows": missing,
            "multiplicity_scope": correction,
            "interpretation": "finite sampled randomization screen conditional on retained data and null mechanism",
        })
    write_tsv(output / "finite_null_resolution.tsv", resolution)

    sensitivity = []
    for threshold in (0.05, 0.10, 0.20):
        spatial = sum(bool(r["spatial_null_p"]) and float(r["spatial_null_p"]) <= threshold for r in score_rows)
        label = sum(bool(r["label_null_p"]) and float(r["label_null_p"]) <= threshold for r in score_rows)
        fake = sum(bool(r["fake_lr_null_p"]) and float(r["fake_lr_null_p"]) <= threshold for r in score_rows)
        all_pass = sum(
            all(r[field] and float(r[field]) <= threshold for field in ["spatial_null_p", "label_null_p", "fake_lr_null_p"])
            for r in score_rows
        )
        high = sum(
            pilot_tier(r, float(r["spatiallr_trust_score_pilot_v3"]), threshold) == "high_pilot_v3"
            for r in score_rows
        )
        sensitivity.append({
            "nominal_threshold": threshold,
            "spatial_attainable_cutoff": math.floor(threshold * 51) / 51,
            "label_attainable_cutoff": math.floor(threshold * 31) / 31,
            "fake_attainable_cutoff": math.floor(threshold * 101) / 101,
            "spatial_pass_rows": spatial,
            "label_pass_rows": label,
            "fake_pass_rows": fake,
            "all_three_pass_rows": all_pass,
            "high_tier_rows_under_gate": high,
            "missing_fake_rows_fail_closed": 9,
        })
    write_tsv(output / "pilot_null_threshold_sensitivity.tsv", sensitivity)

    dataset_samples: dict[str, set[str]] = defaultdict(set)
    passing: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in score_rows:
        dataset_samples[row["dataset"]].add(row["sample_id"])
        if float(row["spatial_null_p"]) <= 0.10 or float(row["label_null_p"]) <= 0.10:
            passing[tuple(row[field] for field in RECURRENCE_KEY)].add(row["sample_id"])
    detail = []
    transitions: Counter[tuple[str, str]] = Counter()
    for row in score_rows:
        key = tuple(row[field] for field in RECURRENCE_KEY)
        contributors = passing.get(key, set())
        if len(contributors) != int(row["dataset_recurrence_n"]):
            raise ValueError("Pilot recurrence reconstruction differs from frozen recurrence count")
        total = len(dataset_samples[row["dataset"]])
        held_out_pass = int(row["sample_id"] in contributors)
        train_pass = len(contributors) - held_out_pass
        train_total = total - 1
        loo_fraction = train_pass / train_total
        loo_support = min(1.0, 2.0 * loo_fraction)
        loo_score = (
            0.40 * float(row["null_support_mean"])
            + 0.15 * float(row["pilot_score_percentile"])
            + 0.15 * loo_support
            + 0.10 * float(row["annotation_quality"])
            + 0.20 * float(row["external_method_support_fraction"])
        )
        loo_tier = pilot_tier(row, loo_score, 0.10)
        original_tier = row["confidence_tier_pilot_v3"]
        transitions[(original_tier, loo_tier)] += 1
        detail.append({
            **{field: row[field] for field in KEY_FIELDS},
            "recurrence_key": "|".join(key),
            "eligible_sample_ids": ";".join(sorted(dataset_samples[row["dataset"]])),
            "biological_unit_kind": "public section/sample_id; no technical replicates in pilot table",
            "passing_sample_ids": ";".join(sorted(contributors)),
            "held_out_sample_passed": held_out_pass,
            "training_samples_total": train_total,
            "training_samples_passed": train_pass,
            "loo_recurrence_fraction": loo_fraction,
            "loo_recurrence_support": loo_support,
            "original_score": row["spatiallr_trust_score_pilot_v3"],
            "loo_score": loo_score,
            "original_tier": original_tier,
            "loo_tier": loo_tier,
            "score_delta": loo_score - float(row["spatiallr_trust_score_pilot_v3"]),
            "tier_changed": int(original_tier != loo_tier),
        })
    write_tsv(output / "pilot_v3_loo_recurrence.tsv", detail)
    summary_rows = [
        {"metric": "rows", "value": len(detail)},
        {"metric": "original_high", "value": sum(r["original_tier"] == "high_pilot_v3" for r in detail)},
        {"metric": "loo_high", "value": sum(r["loo_tier"] == "high_pilot_v3" for r in detail)},
        {"metric": "tier_changed", "value": sum(r["tier_changed"] for r in detail)},
    ]
    for (before, after), count in sorted(transitions.items()):
        summary_rows.append({"metric": f"transition:{before}->{after}", "value": count})
    write_tsv(output / "pilot_v3_loo_summary.tsv", summary_rows)
    return {
        "pilot_rows": len(detail),
        "pilot_original_high": sum(r["original_tier"] == "high_pilot_v3" for r in detail),
        "pilot_loo_high": sum(r["loo_tier"] == "high_pilot_v3" for r in detail),
        "pilot_tier_changed": sum(r["tier_changed"] for r in detail),
    }


def build_f5_comparison_and_loo(project: Path, output: Path) -> dict[str, object]:
    path = project / "scores/spatiallr_trust_score_f5_semisim_full.tsv"
    frame = pd.read_csv(path, sep="\t")
    labelled = frame[frame["resource_axis"].eq(1) & frame["truth_class"].isin(["positive", "negative"])].copy()
    common = labelled[labelled["method_scope_count"].eq(3)].copy()
    if len(frame) != 1140 or len(labelled) != 1020 or len(common) != 510:
        raise ValueError("Unexpected F5 score cardinality")
    units = common.groupby(["sample_id", "replicate"]).size()
    if len(units) != 15 or set(units) != {34}:
        raise ValueError("F5 common-universe sample-replicate structure differs from 15 x 34")
    if set(frame["sample_id"].astype(str)) != {"GSM7780153", "GSM7780154", "GSM7780155"}:
        raise ValueError("F5 biological-section identities differ from the frozen set")
    for _, section in frame.groupby("sample_id"):
        if set(section["replicate"].astype(int)) != {1, 2, 3, 4, 5}:
            raise ValueError("F5 technical replicate set differs from 1-5")
    methods = {
        "COMMOT": ("commot_score", "commot_scoped", "commot_detected"),
        "LIANA": ("liana_score", "liana_scoped", "liana_detected"),
        "CellChat": ("cellchat_score", "cellchat_scoped", "cellchat_detected"),
        "SpatialLR-Trust": ("spatiallr_trust_score_f5_full", "score_eligible", None),
    }
    comparability = []
    for method, (score_col, scoped_col, detected_col) in methods.items():
        scoped = labelled[labelled[scoped_col].eq(1)]
        if detected_col:
            selected = scoped[scoped[detected_col].eq(1)]
            native_rule = f"{detected_col}=1"
        else:
            selected = scoped[scoped["confidence_tier_f5_full"].isin(["high_f5_full", "medium_f5_full"])]
            native_rule = "frozen non-low tier"
        fp = int(selected["truth_class"].eq("negative").sum())
        comparability.append({
            "scorer": method,
            "score_column": score_col,
            "native_callable_rows": len(scoped),
            "common_callable_rows": len(common),
            "native_selection_rule": native_rule,
            "native_selected_rows": len(selected),
            "native_true_positive_rows": int(selected["truth_class"].eq("positive").sum()),
            "native_false_positive_rows": fp,
            "native_empirical_fdp": fp / len(selected) if len(selected) else 0.0,
            "comparison_boundary": "method-native scope is descriptive; head-to-head metrics use only common callable rows",
        })
    write_tsv(output / "f5_method_comparability.tsv", comparability)

    matched = []
    pooled_selected: dict[str, list[pd.DataFrame]] = defaultdict(list)
    for (sample_id, replicate), group in common.groupby(["sample_id", "replicate"], sort=True):
        for method, (score_col, _, _) in methods.items():
            ordered = group.sort_values(
                [score_col, "condition_id", "axis_id"],
                ascending=[False, True, True],
                kind="mergesort",
            )
            selected = ordered.head(4)
            pooled_selected[method].append(selected)
            tp = int(selected["truth_class"].eq("positive").sum())
            matched.append({
                "universe": "common_three_method_labelled",
                "sample_id": sample_id,
                "replicate": int(replicate),
                "scorer": method,
                "score_column": score_col,
                "callable_rows": len(group),
                "top_k": 4,
                "selected_rows": len(selected),
                "true_positive_rows": tp,
                "false_positive_rows": len(selected) - tp,
                "precision_at_k": tp / len(selected),
                "recall_at_k": tp / int(group["truth_class"].eq("positive").sum()),
                "tie_break": "condition_id,axis_id",
                "biological_unit": "sample_id; technical replicate retained within sample",
            })
    pooled_metrics = {}
    labels = common["truth_class"].eq("positive").astype(int)
    for method, (score_col, _, _) in methods.items():
        selected = pd.concat(pooled_selected[method], ignore_index=True)
        tp = int(selected["truth_class"].eq("positive").sum())
        matched.append({
            "universe": "common_three_method_labelled",
            "sample_id": "all",
            "replicate": "all",
            "scorer": method,
            "score_column": score_col,
            "callable_rows": len(common),
            "top_k": "4_per_sample_replicate",
            "selected_rows": len(selected),
            "true_positive_rows": tp,
            "false_positive_rows": len(selected) - tp,
            "precision_at_k": tp / len(selected),
            "recall_at_k": tp / int(labels.sum()),
            "tie_break": "condition_id,axis_id",
            "biological_unit": "3 sections; 5 technical replicates per section",
        })
        pooled_metrics[method] = {
            "auroc": roc_auc(labels, common[score_col]),
            "average_precision": average_precision(labels, common[score_col]),
            "top4_true_positive": tp,
        }
    write_tsv(output / "f5_matched_operating_points.tsv", matched)

    loo_rows = []
    transitions: Counter[tuple[str, str]] = Counter()
    for row in frame.to_dict("records"):
        if int(row["resource_axis"]) != 1:
            loo_rows.append({
                **{field: row[field] for field in F5_ROW_KEY},
                "score_eligible": 0,
                "original_tier": row["confidence_tier_f5_full"],
                "loo_tier": "out_of_scope_fixed_fake",
                "tier_changed": 0,
            })
            transitions[(row["confidence_tier_f5_full"], "out_of_scope_fixed_fake")] += 1
            continue
        passed_ids = set(str(row["passed_biological_sample_ids"]).split(";")) if str(row["passed_biological_sample_ids"]) not in {"", "nan"} else set()
        total = int(row["biological_samples_total"])
        held_out_pass = int(str(row["sample_id"]) in passed_ids)
        train_pass = int(row["joint_pair_null_pass_biological_samples"]) - held_out_pass
        train_total = total - 1
        loo_fraction = train_pass / train_total
        loo_score = (
            0.50 * float(row["pair_null_support_fraction"])
            + 0.20 * float(row["method_support_fraction_all"])
            + 0.20 * loo_fraction
            + 0.10 * float(row["global_null_support_fraction"])
        )
        high = (
            loo_score >= 0.75
            and float(row["observed_local_lr_score"]) > 0
            and int(row["pair_null_pass_count"]) == 3
            and int(row["methods_detected"]) >= 2
            and train_pass >= 2
        )
        medium = (
            not high and loo_score >= 0.50 and float(row["observed_local_lr_score"]) > 0
            and int(row["pair_null_pass_count"]) >= 2 and int(row["methods_detected"]) >= 1
        )
        loo_tier = "high_f5_full" if high else "medium_f5_full" if medium else "low_f5_full"
        original_tier = row["confidence_tier_f5_full"]
        transitions[(original_tier, loo_tier)] += 1
        loo_rows.append({
            **{field: row[field] for field in F5_ROW_KEY},
            "score_eligible": 1,
            "biological_unit_kind": "Xenium section/sample_id; all five technical replicates excluded together",
            "biological_samples_total": total,
            "passed_biological_sample_ids": ";".join(sorted(passed_ids)),
            "held_out_sample_passed": held_out_pass,
            "training_samples_total": train_total,
            "training_samples_passed": train_pass,
            "loo_recurrence_fraction": loo_fraction,
            "original_score": row["spatiallr_trust_score_f5_full"],
            "loo_score": loo_score,
            "score_delta": loo_score - float(row["spatiallr_trust_score_f5_full"]),
            "original_tier": original_tier,
            "loo_tier": loo_tier,
            "tier_changed": int(original_tier != loo_tier),
        })
    write_tsv(output / "f5_loo_recurrence.tsv", loo_rows)
    loo_frame = pd.DataFrame(loo_rows)
    eligible_loo = loo_frame[loo_frame["score_eligible"].eq(1)].copy()
    summary_rows = [
        {"metric": "rows", "value": len(loo_rows)},
        {"metric": "eligible_rows", "value": len(eligible_loo)},
        {"metric": "original_high", "value": int(eligible_loo["original_tier"].eq("high_f5_full").sum())},
        {"metric": "loo_high", "value": int(eligible_loo["loo_tier"].eq("high_f5_full").sum())},
        {"metric": "tier_changed", "value": int(eligible_loo["tier_changed"].sum())},
    ]
    for (before, after), count in sorted(transitions.items()):
        summary_rows.append({"metric": f"transition:{before}->{after}", "value": count})
    metrics = {}
    for name, base in [("common", common), ("all_resource", labelled)]:
        join = base[F5_ROW_KEY + ["spatiallr_trust_score_f5_full", "truth_class"]].merge(
            eligible_loo[F5_ROW_KEY + ["loo_score"]], on=F5_ROW_KEY, validate="one_to_one"
        )
        y = join["truth_class"].eq("positive").astype(int)
        metrics[name] = {
            "original_auroc": roc_auc(y, join["spatiallr_trust_score_f5_full"]),
            "loo_auroc": roc_auc(y, join["loo_score"]),
            "original_average_precision": average_precision(y, join["spatiallr_trust_score_f5_full"]),
            "loo_average_precision": average_precision(y, join["loo_score"]),
        }
        for metric, value in metrics[name].items():
            summary_rows.append({"metric": f"{name}:{metric}", "value": value})
    write_tsv(output / "f5_loo_summary.tsv", summary_rows)
    return {
        "f5_rows": len(frame),
        "f5_common_rows": len(common),
        "f5_original_high": int(eligible_loo["original_tier"].eq("high_f5_full").sum()),
        "f5_loo_high": int(eligible_loo["loo_tier"].eq("high_f5_full").sum()),
        "f5_tier_changed": int(eligible_loo["tier_changed"].sum()),
        "matched_metrics": pooled_metrics,
        "loo_metrics": metrics,
    }


def open_text(path: Path):
    return gzip.open(path, "rt", errors="replace") if path.suffix == ".gz" else path.open("rt", errors="replace")


def target_metrics(
    genes: list[str], receiver: str, comp_counts: Counter[str],
    gene_sum: dict[str, Counter[str]], gene_spots: dict[str, dict[str, set[int]]],
) -> dict[str, object]:
    detected = [gene for gene in genes if gene in gene_sum]
    detected_n = len(detected)
    comp_mean = {}
    comp_positive = {}
    for comp in COMPARTMENTS:
        spots = comp_counts.get(comp, 0)
        total = sum(gene_sum[gene][comp] for gene in detected)
        comp_mean[comp] = total / max(spots * max(detected_n, 1), 1)
        positive: set[int] = set()
        for gene in detected:
            positive.update(gene_spots.get(gene, {}).get(comp, set()))
        comp_positive[comp] = len(positive) / spots if spots else 0.0
    other_spots = sum(comp_counts.get(comp, 0) for comp in COMPARTMENTS if comp != receiver)
    other_total = sum(gene_sum[gene][comp] for gene in detected for comp in COMPARTMENTS if comp != receiver)
    other_mean = other_total / max(other_spots * max(detected_n, 1), 1)
    receiver_mean = comp_mean.get(receiver, 0.0)
    percentile = sum(value <= receiver_mean for value in sorted(comp_mean.values())) / len(COMPARTMENTS)
    enrichment = (receiver_mean + 1e-9) / (other_mean + 1e-9)
    support = int(detected_n >= 2 and percentile >= 0.75 and enrichment >= 1.25)
    return {
        "target_genes_in_panel": len(genes),
        "target_genes_detected": detected_n,
        "target_genes_detected_list": ",".join(detected),
        "receiver_target_mean": receiver_mean,
        "other_compartment_target_mean": other_mean,
        "receiver_target_enrichment": enrichment,
        "receiver_target_positive_spot_fraction": comp_positive.get(receiver, 0.0),
        "receiver_target_percentile": percentile,
        "target_support": support,
    }


def load_sample_targets(
    project: Path, row: dict[str, str], annotation: dict[str, str], wanted: set[str], scratch: Path,
) -> tuple[Counter[str], dict[str, Counter[str]], dict[str, dict[str, set[int]]]]:
    sample_dir = scratch / f"{row['dataset']}__{row['sample_id']}"
    sample_dir.mkdir(parents=True, exist_ok=False)
    local = {}
    for role, source in [
        ("barcodes", project / row["barcodes_path"]),
        ("features", project / row["features_path"]),
        ("expression", project / row["expression_path"]),
        ("annotation", project / annotation["annotation_path"]),
    ]:
        destination = sample_dir / source.name
        shutil.copy2(source, destination)
        local[role] = destination
    try:
        with open_text(local["barcodes"]) as handle:
            barcodes = [line.rstrip("\n").split("\t", 1)[0] for line in handle if line.strip()]
        with gzip.open(local["annotation"], "rt", errors="replace", newline="") as handle:
            barcode_to_compartment = {r["barcode"]: r["spot_compartment"] for r in csv.DictReader(handle, delimiter="\t")}
        comp_by_index = {index: barcode_to_compartment.get(barcode) for index, barcode in enumerate(barcodes)}
        comp_counts = Counter(comp for comp in comp_by_index.values() if comp in COMPARTMENTS)
        gene_index = {}
        with open_text(local["features"]) as handle:
            for index, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                names = [value.upper() for value in fields[:2]]
                hit = next((name for name in names if name in wanted), None)
                if hit:
                    gene_index[index] = hit
        gene_sum: dict[str, Counter[str]] = defaultdict(Counter)
        gene_spots: dict[str, dict[str, set[int]]] = defaultdict(lambda: defaultdict(set))
        shape_seen = False
        with open_text(local["expression"]) as handle:
            for line in handle:
                if line.startswith("%"):
                    continue
                fields = line.split()
                if not shape_seen:
                    shape_seen = True
                    continue
                gene = gene_index.get(int(fields[0]))
                comp = comp_by_index.get(int(fields[1]) - 1)
                if gene and comp in COMPARTMENTS:
                    value = float(fields[2])
                    gene_sum[gene][comp] += value
                    gene_spots[gene][comp].add(int(fields[1]) - 1)
        return comp_counts, dict(gene_sum), {gene: dict(value) for gene, value in gene_spots.items()}
    finally:
        shutil.rmtree(sample_dir)


def build_target_sensitivity(project: Path, output: Path, scratch: Path) -> dict[str, int]:
    panel_rows = read_tsv(project / "docs/pilot_pathway_target_panel.tsv")
    pathway_genes: dict[str, list[str]] = defaultdict(list)
    for row in panel_rows:
        gene = row["gene"].upper()
        if gene not in pathway_genes[row["pathway"]]:
            pathway_genes[row["pathway"]].append(gene)
    wanted = {gene for genes in pathway_genes.values() for gene in genes}
    score_rows = read_tsv(project / "scores/spatiallr_trust_score_pilot_v3.tsv")
    score_keys = [tuple(row[field] for field in KEY_FIELDS) for row in score_rows]
    if len(score_keys) != len(set(score_keys)):
        raise ValueError("Pilot score keys are not unique")
    current_source = read_tsv(project / "results/task_e_target_activation/target_activation_all_score_v3_with_target.tsv")
    current_keys = [tuple(row[field] for field in KEY_FIELDS) for row in current_source]
    if len(current_keys) != len(set(current_keys)) or set(current_keys) != set(score_keys):
        raise ValueError("Current target table is not an exact one-to-one pilot-score join")
    current_rows = {tuple(row[field] for field in KEY_FIELDS): row for row in current_source}
    candidates_by_sample: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in score_rows:
        candidates_by_sample[(row["dataset"], row["sample_id"])].append(row)
    manifests = [row for row in read_tsv(project / "docs/minimal_input_manifest.tsv") if row.get("status", "ready") == "ready"]
    manifest_keys = [(row["dataset"], row["sample_id"]) for row in manifests]
    if len(manifest_keys) != len(set(manifest_keys)) or len(manifests) != 34:
        raise ValueError("Ready minimal-input manifest must contain 34 unique dataset/sample rows")
    annotation_source = read_tsv(project / "docs/spot_annotation_manifest.tsv")
    annotation_keys = [(row["dataset"], row["sample_id"]) for row in annotation_source]
    if len(annotation_keys) != len(set(annotation_keys)):
        raise ValueError("Spot-annotation manifest keys are not unique")
    annotations = {(row["dataset"], row["sample_id"]): row for row in annotation_source}
    if not set(manifest_keys).issubset(annotations):
        raise ValueError("Ready minimal inputs do not have one annotation row each")

    candidate_genes: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    candidate_pairs: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in score_rows:
        for role in ("ligand", "receptor"):
            gene = row[role].upper()
            if gene in pathway_genes[row["pathway"]]:
                candidate_genes[row["pathway"]][gene].add(role)
                candidate_pairs[(row["pathway"], gene)].add(f"{row['ligand'].upper()}-{row['receptor'].upper()}")
    overlap = []
    for pathway in sorted(candidate_genes):
        for gene in sorted(candidate_genes[pathway]):
            overlap.append({
                "pathway": pathway,
                "overlap_gene": gene,
                "candidate_roles": ";".join(sorted(candidate_genes[pathway][gene])),
                "candidate_pairs": ";".join(sorted(candidate_pairs[(pathway, gene)])),
                "interpretation": "candidate-defining gene overlaps the pathway co-expression panel",
            })
    write_tsv(output / "target_panel_overlap.tsv", overlap)

    detail = []
    mif_loo = []
    reproduced_support = 0
    reproduced_high_support = 0
    seen_mif = set()
    sufficient_statistics = []
    float_fields = [
        "receiver_target_mean", "other_compartment_target_mean", "receiver_target_enrichment",
        "receiver_target_positive_spot_fraction", "receiver_target_percentile",
    ]
    for manifest in manifests:
        sample_id = manifest["sample_id"]
        comp_counts, gene_sum, gene_spots = load_sample_targets(
            project, manifest, annotations[(manifest["dataset"], sample_id)], wanted, scratch
        )
        for gene in sorted(wanted):
            for comp in COMPARTMENTS:
                positive_indices = gene_spots.get(gene, {}).get(comp, set())
                positive_bitset = 0
                for index in positive_indices:
                    positive_bitset |= 1 << index
                sufficient_statistics.append({
                    "dataset": manifest["dataset"],
                    "sample_id": sample_id,
                    "gene": gene,
                    "compartment": comp,
                    "compartment_spots": comp_counts.get(comp, 0),
                    "gene_expression_sum": gene_sum.get(gene, Counter()).get(comp, 0.0),
                    "gene_positive_spots": len(positive_indices),
                    "positive_spot_bitset_hex": format(positive_bitset, "x"),
                    "gene_detected_in_matrix": int(gene in gene_sum),
                })
        for row in candidates_by_sample[(manifest["dataset"], sample_id)]:
            genes = pathway_genes[row["pathway"]]
            original = target_metrics(genes, row["receiver_compartment"], comp_counts, gene_sum, gene_spots)
            removed = [gene for gene in genes if gene in {row["ligand"].upper(), row["receptor"].upper()}]
            filtered_genes = [gene for gene in genes if gene not in set(removed)]
            excluded = target_metrics(filtered_genes, row["receiver_compartment"], comp_counts, gene_sum, gene_spots)
            current = current_rows[tuple(row[field] for field in KEY_FIELDS)]
            if int(current["target_activation_support"] or 0) != original["target_support"]:
                raise ValueError("Original target support did not reproduce")
            for field in float_fields:
                if not math.isclose(float(current[field]), float(original[field]), rel_tol=2e-6, abs_tol=2e-8):
                    raise ValueError(f"Original target metric did not reproduce: {field}")
            reproduced_support += int(original["target_support"])
            if row["confidence_tier_pilot_v3"] == "high_pilot_v3":
                reproduced_high_support += int(original["target_support"])
            out = {field: row[field] for field in KEY_FIELDS}
            out.update({
                "confidence_tier_pilot_v3": row["confidence_tier_pilot_v3"],
                "candidate_genes_removed": ",".join(removed),
            })
            out.update({f"original_{key}": value for key, value in original.items()})
            out.update({f"excluded_{key}": value for key, value in excluded.items()})
            out.update({
                "support_flip": f"{original['target_support']}->{excluded['target_support']}",
                "enrichment_delta": float(excluded["receiver_target_enrichment"]) - float(original["receiver_target_enrichment"]),
                "interpretation": "candidate-excluded receiver-compartment target-gene co-expression proxy",
            })
            detail.append(out)
            if row["pathway"] == "MIF":
                identity = (row["dataset"], sample_id, row["receiver_compartment"])
                if identity not in seen_mif:
                    seen_mif.add(identity)
                    for omitted in genes:
                        metric = target_metrics([gene for gene in genes if gene != omitted], row["receiver_compartment"], comp_counts, gene_sum, gene_spots)
                        mif_loo.append({
                            "dataset": row["dataset"], "sample_id": sample_id,
                            "receiver_compartment": row["receiver_compartment"], "omitted_gene": omitted,
                            **metric,
                            "featured_case": int(sample_id == "GSM8855712_NSCLC_P10" and row["receiver_compartment"] == "tumor_like"),
                        })
    if len(detail) != 4225 or reproduced_support != 1266 or reproduced_high_support != 134:
        raise ValueError("Target proxy reproduction cardinality mismatch")
    write_tsv(output / "target_candidate_excluded.tsv", detail)
    write_tsv(output / "mif_target_gene_loo.tsv", mif_loo)
    write_tsv(output / "target_gene_compartment_sufficient_statistics.tsv", sufficient_statistics)
    summary_rows = [
        {"metric": "rows", "value": len(detail)},
        {"metric": "overlap_genes", "value": len(overlap)},
        {"metric": "original_supported_rows", "value": reproduced_support},
        {"metric": "excluded_supported_rows", "value": sum(int(row["excluded_target_support"]) for row in detail)},
        {"metric": "support_flips", "value": sum(row["original_target_support"] != row["excluded_target_support"] for row in detail)},
        {"metric": "high_original_supported_rows", "value": reproduced_high_support},
        {"metric": "high_excluded_supported_rows", "value": sum(
            int(row["excluded_target_support"]) for row in detail if row["confidence_tier_pilot_v3"] == "high_pilot_v3"
        )},
        {"metric": "mif_rows", "value": sum(row["pathway"] == "MIF" for row in detail)},
        {"metric": "mif_original_supported_rows", "value": sum(
            int(row["original_target_support"]) for row in detail if row["pathway"] == "MIF"
        )},
        {"metric": "mif_excluded_supported_rows", "value": sum(
            int(row["excluded_target_support"]) for row in detail if row["pathway"] == "MIF"
        )},
    ]
    write_tsv(output / "target_candidate_excluded_summary.tsv", summary_rows)
    return {row["metric"]: int(row["value"]) for row in summary_rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--scratch-dir", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    project = Path(args.project).resolve()
    output = Path(args.output_dir).resolve()
    scratch = Path(args.scratch_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=True)
    failures: list[dict[str, str]] = []
    try:
        construct_rows = build_construct_provenance(project, output, args.source_commit)
        pilot = build_null_and_pilot_loo(project, output)
        f5 = build_f5_comparison_and_loo(project, output)
        target = build_target_sensitivity(project, output, scratch)
        summary = {
            "status": "PASS_TF48_STAGE4_REVISION_ANALYSES_BUILT",
            "run_id": args.run_id,
            "source_commit": args.source_commit,
            "construct_rows": construct_rows,
            **pilot,
            **f5,
            "target": target,
            "frozen_scores_modified": False,
            "interpretation": "Stage 4 computational sensitivity and provenance package; not biological validation or calibrated probability.",
            "validation_pending": True,
        }
        (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    except Exception as exc:
        failures.append({"error": repr(exc)})
        raise
    finally:
        write_tsv(output / "failures.tsv", failures, ["error"])


if __name__ == "__main__":
    main()
