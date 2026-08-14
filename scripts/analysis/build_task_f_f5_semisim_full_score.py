#!/usr/bin/env python3
"""Build the frozen F5 full semi-simulation score and descriptive metrics."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROW_KEY = [
    "sample_id",
    "replicate",
    "condition_id",
    "axis_id",
    "sender_compartment",
    "receiver_compartment",
    "ligand",
    "receptor",
]
RECURRENCE_KEY = [
    "condition_id",
    "sender_compartment",
    "receiver_compartment",
    "ligand",
    "receptor",
]
PAIR_Q = [
    "coordinate_pair_bh_fdr",
    "label_pair_bh_fdr",
    "fake_lr_pair_bh_fdr",
]
GLOBAL_Q = ["coordinate_bh_fdr", "label_bh_fdr", "fake_lr_bh_fdr"]
METHODS = {
    "commot": ("commot_score", "commot_output", "commot_detected"),
    "liana": ("liana_score", "liana_output", "liana_detected"),
    "cellchat": ("cellchat_score", "cellchat_output", "cellchat_detected"),
}
RESOURCE_NAMES = {
    "commot": "COMMOT",
    "liana": "LIANA",
    "cellchat": "CellChat",
}
SCORE = "spatiallr_trust_score_f5_full"
TIER = "confidence_tier_f5_full"
TOP_K = 4
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20260726


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require_unique(frame: pd.DataFrame, keys: list[str], name: str) -> None:
    if frame.duplicated(keys).any():
        raise ValueError(f"{name} contains duplicate keys: {keys}")


def condition_family(condition_id: str) -> str:
    if condition_id.startswith("coordinate_mix_"):
        return "coordinate_corruption"
    if condition_id.startswith("label_mix_"):
        return "label_corruption"
    if condition_id.startswith("expression_mix_"):
        return "expression_corruption"
    if condition_id.startswith("spike_"):
        return "spike_strength"
    if condition_id == "unmodified_reference":
        return "reference"
    if condition_id.startswith("hard_negative_"):
        return "hard_negative"
    if condition_id == "fixed_fake_lr":
        return "fixed_fake_control"
    raise ValueError(f"Unknown condition: {condition_id}")


def roc_auc(labels: pd.Series, scores: pd.Series) -> float:
    labels_np = labels.astype(int).to_numpy()
    positives = labels_np == 1
    n_pos = int(positives.sum())
    n_neg = int((~positives).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = scores.astype(float).rank(method="average").to_numpy()
    return float(
        (ranks[positives].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    )


def average_precision(labels: pd.Series, scores: pd.Series) -> float:
    frame = pd.DataFrame(
        {"label": labels.astype(int), "score": scores.astype(float)}
    )
    if frame["label"].sum() == 0:
        return float("nan")
    grouped = (
        frame.groupby("score", sort=False)["label"]
        .agg(["sum", "count"])
        .sort_index(ascending=False)
    )
    true_positive = grouped["sum"].cumsum()
    false_positive = (grouped["count"] - grouped["sum"]).cumsum()
    recall = true_positive / frame["label"].sum()
    precision = true_positive / (true_positive + false_positive)
    return float(((recall - recall.shift(fill_value=0)) * precision).sum())


def metric_row(
    frame: pd.DataFrame,
    universe: str,
    scorer: str,
    score_column: str,
    sample_id: str = "all",
) -> dict[str, object]:
    ordered = frame.sort_values(
        [score_column, "sample_id", "replicate", "condition_id", "axis_id"],
        ascending=[False, True, True, True, True],
        kind="mergesort",
    )
    labels = frame["truth_class"].eq("positive").astype(int)
    top = ordered.head(TOP_K)
    positives = int(labels.sum())
    selected = len(top)
    top_positive = int(top["truth_class"].eq("positive").sum())
    return {
        "sample_id": sample_id,
        "metric_universe": universe,
        "scorer": scorer,
        "score_column": score_column,
        "rows": len(frame),
        "positives": positives,
        "negatives": int((1 - labels).sum()),
        "auroc": roc_auc(labels, frame[score_column]),
        "average_precision": average_precision(labels, frame[score_column]),
        "top_k": TOP_K,
        "top_k_selected": selected,
        "top_k_positives": top_positive,
        "precision_at_k": top_positive / selected if selected else float("nan"),
        "recall_at_k": top_positive / positives if positives else float("nan"),
        "tie_break": "sample_replicate_condition_axis",
    }


def metric_specs(score: pd.DataFrame) -> list[tuple[str, str, str, pd.DataFrame]]:
    labelled = score[
        score["resource_axis"].eq(1)
        & score["truth_class"].isin(["positive", "negative"])
    ].copy()
    common = labelled[labelled["method_scope_count"].eq(3)]
    specs: list[tuple[str, str, str, pd.DataFrame]] = []
    for method, (column, _, _) in METHODS.items():
        specs.append(("common_three_method_primary", method, column, common))
    specs.append(("common_three_method_primary", "spatiallr_trust", SCORE, common))
    for method, (column, output, _) in METHODS.items():
        specs.append(
            (
                "method_scoped_resource",
                method,
                column,
                labelled[labelled[f"{method}_scoped"].eq(1)],
            )
        )
    specs.append(("all_resource_primary", "spatiallr_trust", SCORE, labelled))
    return specs


def build_metrics(
    score: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pooled_rows = []
    sample_rows = []
    bootstrap_rows = []
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    for universe, scorer, column, frame in metric_specs(score):
        pooled_rows.append(metric_row(frame, universe, scorer, column))
        for sample_id, section in frame.groupby("sample_id", sort=True):
            sample_rows.append(
                metric_row(section, universe, scorer, column, str(sample_id))
            )
        sample_ids = sorted(frame["sample_id"].unique())
        boot_auc = []
        boot_ap = []
        for _ in range(BOOTSTRAP_DRAWS):
            selected_ids = rng.choice(sample_ids, size=len(sample_ids), replace=True)
            pieces = []
            for copy_id, sample_id in enumerate(selected_ids):
                piece = frame[frame["sample_id"].eq(sample_id)].copy()
                piece["_bootstrap_copy"] = copy_id
                pieces.append(piece)
            boot = pd.concat(pieces, ignore_index=True)
            labels = boot["truth_class"].eq("positive").astype(int)
            boot_auc.append(roc_auc(labels, boot[column]))
            boot_ap.append(average_precision(labels, boot[column]))
        for metric_name, values in (
            ("auroc", boot_auc),
            ("average_precision", boot_ap),
        ):
            values_np = np.asarray(values, dtype=float)
            valid = values_np[np.isfinite(values_np)]
            bootstrap_rows.append(
                {
                    "metric_universe": universe,
                    "scorer": scorer,
                    "metric": metric_name,
                    "point_estimate": (
                        pooled_rows[-1]["auroc"]
                        if metric_name == "auroc"
                        else pooled_rows[-1]["average_precision"]
                    ),
                    "ci_lower_0_025": float(np.quantile(valid, 0.025)),
                    "ci_upper_0_975": float(np.quantile(valid, 0.975)),
                    "bootstrap_draws": BOOTSTRAP_DRAWS,
                    "valid_draws": len(valid),
                    "cluster_unit": "sample_id",
                    "biological_samples": len(sample_ids),
                    "seed": BOOTSTRAP_SEED,
                }
            )
    return (
        pd.DataFrame(pooled_rows),
        pd.DataFrame(sample_rows),
        pd.DataFrame(bootstrap_rows),
    )


def build_topk(score: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for universe, scorer, column, frame in metric_specs(score):
        for (sample_id, replicate), part in frame.groupby(
            ["sample_id", "replicate"], sort=True
        ):
            row = metric_row(
                part, universe, scorer, column, sample_id=str(sample_id)
            )
            row["replicate"] = int(replicate)
            rows.append(row)
    return pd.DataFrame(rows)


def build_selection_audits(
    score: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    resource = score[score["resource_axis"].eq(1)].copy()
    labelled = resource[
        resource["truth_class"].isin(["positive", "negative"])
    ].copy()
    selections = {
        "commot_detected": resource["commot_detected"].eq(1),
        "liana_detected": resource["liana_detected"].eq(1),
        "cellchat_detected": resource["cellchat_detected"].eq(1),
        "trust_nonlow": resource[TIER].isin(
            ["high_f5_full", "medium_f5_full"]
        ),
        "trust_high": resource[TIER].eq("high_f5_full"),
    }
    fdp_rows = []
    family_rows = []
    condition_rows = []
    for selection_name, selected_mask in selections.items():
        selected = labelled[selected_mask.loc[labelled.index]]
        false_positive = int(selected["truth_class"].eq("negative").sum())
        fdp_rows.append(
            {
                "selection": selection_name,
                "selected": len(selected),
                "true_positive": int(
                    selected["truth_class"].eq("positive").sum()
                ),
                "false_positive": false_positive,
                "empirical_false_discovery_proportion": (
                    false_positive / len(selected) if len(selected) else 0.0
                ),
                "fixed_fake_excluded": 1,
                "reference_excluded": 1,
            }
        )
        hard = resource[resource["condition_id"].str.startswith("hard_negative_")]
        hard_selected = selected_mask.loc[hard.index]
        for condition_id, family in hard.groupby("condition_id", sort=True):
            index = family.index
            count = int(selected_mask.loc[index].sum())
            family_rows.append(
                {
                    "condition_id": condition_id,
                    "selection": selection_name,
                    "negative_rows": len(family),
                    "selected_rows": count,
                    "false_positive_rate": count / len(family),
                    "samples": family["sample_id"].nunique(),
                    "technical_replicates": family[
                        ["sample_id", "replicate"]
                    ].drop_duplicates().shape[0],
                }
            )
        for condition_id, condition in resource.groupby("condition_id", sort=True):
            index = condition.index
            count = int(selected_mask.loc[index].sum())
            truth_classes = sorted(condition["truth_class"].unique())
            if len(truth_classes) != 1:
                raise ValueError(f"Mixed truth labels in {condition_id}")
            condition_rows.append(
                {
                    "condition_id": condition_id,
                    "condition_family": condition_family(condition_id),
                    "truth_class": truth_classes[0],
                    "selection": selection_name,
                    "rows": len(condition),
                    "selected_rows": count,
                    "selection_rate": count / len(condition),
                }
            )
    return (
        pd.DataFrame(fdp_rows),
        pd.DataFrame(family_rows),
        pd.DataFrame(condition_rows),
    )


def build_monotonicity(score: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    series = {
        "spike_strength": [
            "spike_low",
            "spike_medium",
            "spike_high",
        ],
        "coordinate_corruption": [
            "unmodified_reference",
            "coordinate_mix_10",
            "coordinate_mix_25",
            "coordinate_mix_50",
            "coordinate_mix_100",
        ],
        "label_corruption": [
            "unmodified_reference",
            "label_mix_10",
            "label_mix_25",
            "label_mix_50",
        ],
        "expression_corruption": [
            "unmodified_reference",
            "expression_mix_10",
            "expression_mix_25",
            "expression_mix_50",
        ],
    }
    scorers = {
        "commot": ("commot_score", "commot_scoped"),
        "liana": ("liana_score", "liana_scoped"),
        "cellchat": ("cellchat_score", "cellchat_scoped"),
        "spatiallr_trust": (SCORE, "score_eligible"),
    }
    rows = []
    for family, order in series.items():
        subset = score[score["condition_id"].isin(order)].copy()
        for scorer, (column, eligibility) in scorers.items():
            for group_key, part in subset.groupby(
                [
                    "sample_id",
                    "replicate",
                    "axis_id",
                    "sender_compartment",
                    "receiver_compartment",
                    "ligand",
                    "receptor",
                ],
                sort=True,
            ):
                indexed = part.set_index("condition_id")
                if not set(order).issubset(indexed.index):
                    continue
                selected = indexed.loc[order]
                if not selected[eligibility].eq(1).all():
                    continue
                values = selected[column].astype(float).to_numpy()
                differences = np.diff(values)
                expected = differences >= 0 if family == "spike_strength" else differences <= 0
                rows.append(
                    {
                        "perturbation_family": family,
                        "expected_direction": (
                            "nondecreasing"
                            if family == "spike_strength"
                            else "nonincreasing"
                        ),
                        "scorer": scorer,
                        "sample_id": group_key[0],
                        "replicate": group_key[1],
                        "axis_id": group_key[2],
                        "steps": len(differences),
                        "expected_steps": int(expected.sum()),
                        "fully_monotonic": int(expected.all()),
                        "first_value": values[0],
                        "last_value": values[-1],
                    }
                )
    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby(
            ["perturbation_family", "expected_direction", "scorer"],
            as_index=False,
        )
        .agg(
            series=("fully_monotonic", "size"),
            fully_monotonic=("fully_monotonic", "sum"),
            fully_monotonic_fraction=("fully_monotonic", "mean"),
            expected_step_fraction=(
                "expected_steps",
                lambda values: values.sum()
                / detail.loc[values.index, "steps"].sum(),
            ),
        )
    )
    return detail, summary


def build_rank_stability(score: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scorers = {
        "commot": ("commot_score", "commot_scoped"),
        "liana": ("liana_score", "liana_scoped"),
        "cellchat": ("cellchat_score", "cellchat_scoped"),
        "spatiallr_trust": (SCORE, "score_eligible"),
    }
    key = ["condition_id", "axis_id"]
    rows = []
    for sample_id, section in score[score["resource_axis"].eq(1)].groupby(
        "sample_id", sort=True
    ):
        for scorer, (column, eligibility) in scorers.items():
            for left_rep, right_rep in itertools.combinations(range(1, 6), 2):
                left = section[
                    section["replicate"].eq(left_rep)
                    & section[eligibility].eq(1)
                ][key + [column]]
                right = section[
                    section["replicate"].eq(right_rep)
                    & section[eligibility].eq(1)
                ][key + [column]]
                merged = left.merge(
                    right,
                    on=key,
                    how="inner",
                    validate="one_to_one",
                    suffixes=("_left", "_right"),
                )
                rho = merged[f"{column}_left"].corr(
                    merged[f"{column}_right"], method="spearman"
                )
                rows.append(
                    {
                        "sample_id": sample_id,
                        "scorer": scorer,
                        "replicate_left": left_rep,
                        "replicate_right": right_rep,
                        "common_rows": len(merged),
                        "spearman": rho,
                    }
                )
    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby(["sample_id", "scorer"], as_index=False)
        .agg(
            replicate_pairs=("spearman", "size"),
            mean_spearman=("spearman", "mean"),
            median_spearman=("spearman", "median"),
            min_spearman=("spearman", "min"),
            max_spearman=("spearman", "max"),
        )
    )
    return detail, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method-audit", required=True)
    parser.add_argument("--method-validation-summary", required=True)
    parser.add_argument("--null-audit", required=True)
    parser.add_argument("--null-validation-summary", required=True)
    parser.add_argument("--null-root", required=True)
    parser.add_argument("--tf6-development-score", required=True)
    parser.add_argument("--v3-checksums", required=True)
    parser.add_argument("--output-score", required=True)
    parser.add_argument("--output-meta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    method_path = Path(args.method_audit).resolve()
    null_audit_path = Path(args.null_audit).resolve()
    null_root = Path(args.null_root).resolve()
    output = Path(args.output_dir).resolve()
    method_summary = json.loads(Path(args.method_validation_summary).read_text())
    null_summary = json.loads(Path(args.null_validation_summary).read_text())
    if method_summary["status"] != "PASS_F5_SEMISIM_FULL_METHODS":
        raise ValueError("Full method gate has not passed")
    if null_summary["status"] != "PASS_F5_SEMISIM_FULL_NULL999":
        raise ValueError("Full null gate has not passed")

    method = pd.read_csv(method_path, sep="\t")
    null_audit = pd.read_csv(null_audit_path, sep="\t")
    if len(method) != 1140 or len(null_audit) != 1140:
        raise ValueError("Expected 1,140 method and null truth-axis rows")
    require_unique(method, ROW_KEY, "method audit")
    require_unique(null_audit, ROW_KEY, "null audit")
    identity = ROW_KEY + ["truth_class", "resource_axis"]
    if not method[identity].sort_values(ROW_KEY).reset_index(drop=True).equals(
        null_audit[identity].sort_values(ROW_KEY).reset_index(drop=True)
    ):
        raise ValueError("Method and null truth contracts differ")

    method_columns = ROW_KEY + ["truth_class", "resource_axis"]
    for _, columns in METHODS.items():
        method_columns.extend(columns)
    score = method[method_columns].merge(
        null_audit[
            ROW_KEY
            + ["null_output", *PAIR_Q, *GLOBAL_Q]
        ],
        on=ROW_KEY,
        how="left",
        validate="one_to_one",
    )
    null_frames = []
    for (sample_id, replicate, condition_id), part in score[
        score["resource_axis"].eq(1)
    ].groupby(["sample_id", "replicate", "condition_id"], sort=True):
        path = (
            null_root
            / str(sample_id)
            / f"replicate_{int(replicate)}"
            / str(condition_id)
            / "null_candidate_scores.tsv"
        )
        frame = pd.read_csv(path, sep="\t")
        frame["sample_id"] = sample_id
        frame["replicate"] = int(replicate)
        frame["condition_id"] = condition_id
        wanted = part[
            [
                "sender_compartment",
                "receiver_compartment",
                "ligand",
                "receptor",
            ]
        ].drop_duplicates()
        frame = frame.merge(
            wanted,
            on=[
                "sender_compartment",
                "receiver_compartment",
                "ligand",
                "receptor",
            ],
            how="inner",
            validate="many_to_one",
        )
        null_frames.append(
            frame[
                [
                    "sample_id",
                    "replicate",
                    "condition_id",
                    "sender_compartment",
                    "receiver_compartment",
                    "ligand",
                    "receptor",
                    "resource_methods",
                    "method_scope_count",
                    "observed_local_lr_score",
                    *PAIR_Q,
                    *GLOBAL_Q,
                ]
            ]
        )
    raw_null = pd.concat(null_frames, ignore_index=True)
    raw_key = [
        "sample_id",
        "replicate",
        "condition_id",
        "sender_compartment",
        "receiver_compartment",
        "ligand",
        "receptor",
    ]
    require_unique(raw_null, raw_key, "raw truth-axis null rows")
    score = score.merge(
        raw_null,
        on=raw_key,
        how="left",
        validate="many_to_one",
        suffixes=("_audit", ""),
    )
    resource = score["resource_axis"].eq(1)
    if score.loc[resource, "observed_local_lr_score"].isna().any():
        raise ValueError("Resource truth axis missing raw null evidence")
    for column in [*PAIR_Q, *GLOBAL_Q]:
        if not np.allclose(
            score.loc[resource, f"{column}_audit"],
            score.loc[resource, column],
            rtol=0,
            atol=1e-12,
        ):
            raise ValueError(f"Raw/audit mismatch: {column}")
        score[column] = score[column].fillna(score[f"{column}_audit"])
        score = score.drop(columns=f"{column}_audit")
    score["resource_methods"] = score["resource_methods"].fillna("")
    score["method_scope_count"] = score["method_scope_count"].fillna(0).astype(int)
    score["observed_local_lr_score"] = score["observed_local_lr_score"].fillna(0.0)
    score["condition_family"] = score["condition_id"].map(condition_family)

    detected_columns = [columns[2] for columns in METHODS.values()]
    for method, resource_name in RESOURCE_NAMES.items():
        score[f"{method}_scoped"] = score["resource_methods"].map(
            lambda value, name=resource_name: int(name in value.split(";"))
        )
        if (score[METHODS[method][1]] > score[f"{method}_scoped"]).any():
            raise ValueError(f"{method} output exceeds resource scope")
    if not np.array_equal(
        score[[f"{method}_scoped" for method in METHODS]].sum(axis=1).astype(int),
        score["method_scope_count"],
    ):
        raise ValueError("Derived method scope does not match null resource scope")
    score["methods_detected"] = score[detected_columns].sum(axis=1).astype(int)
    if (score["methods_detected"] > score["method_scope_count"]).any():
        raise ValueError("Method detections exceed scope")
    score["method_support_fraction_all"] = score["methods_detected"] / 3
    score["method_support_fraction_scoped"] = np.divide(
        score["methods_detected"],
        score["method_scope_count"],
        out=np.zeros(len(score), dtype=float),
        where=score["method_scope_count"].to_numpy() > 0,
    )
    pair_pass = []
    global_pass = []
    for column in PAIR_Q:
        pass_column = column.replace("_bh_fdr", "_pass")
        score[pass_column] = score[column].le(0.05).astype(np.int8)
        pair_pass.append(pass_column)
    for column in GLOBAL_Q:
        pass_column = column.replace("_bh_fdr", "_pass")
        score[pass_column] = score[column].le(0.05).astype(np.int8)
        global_pass.append(pass_column)
    score["pair_null_pass_count"] = score[pair_pass].sum(axis=1)
    score["pair_null_support_fraction"] = score["pair_null_pass_count"] / 3
    score["global_null_pass_count"] = score[global_pass].sum(axis=1)
    score["global_null_support_fraction"] = score["global_null_pass_count"] / 3
    score["joint_pair_null_pass"] = score["pair_null_pass_count"].eq(3).astype(np.int8)

    technical_key = ["sample_id", *RECURRENCE_KEY]
    technical = (
        score.groupby(technical_key, as_index=False)
        .agg(
            technical_replicates_total=("replicate", "nunique"),
            technical_replicates_joint_pass=("joint_pair_null_pass", "sum"),
            resource_axis=("resource_axis", "min"),
        )
    )
    if set(technical["technical_replicates_total"]) != {5}:
        raise ValueError("Every biological sample-axis must have five replicates")
    technical["biological_sample_joint_pass_all_replicates"] = (
        technical["technical_replicates_joint_pass"]
        .eq(technical["technical_replicates_total"])
        .astype(np.int8)
    )
    recurrence = (
        technical.groupby(RECURRENCE_KEY, as_index=False)
        .agg(
            biological_samples_total=("sample_id", "nunique"),
            biological_samples_joint_pass_all_replicates=(
                "biological_sample_joint_pass_all_replicates",
                "sum",
            ),
            technical_replicates_total=(
                "technical_replicates_total",
                "sum",
            ),
            technical_replicates_joint_pass=(
                "technical_replicates_joint_pass",
                "sum",
            ),
            resource_axis=("resource_axis", "min"),
        )
    )
    if set(recurrence["biological_samples_total"]) != {3}:
        raise ValueError("Recurrence denominator must be three biological samples")
    recurrence["joint_pair_null_recurrence_fraction"] = (
        recurrence["biological_samples_joint_pass_all_replicates"]
        / recurrence["biological_samples_total"]
    )
    passed_ids = (
        technical[technical["biological_sample_joint_pass_all_replicates"].eq(1)]
        .groupby(RECURRENCE_KEY)["sample_id"]
        .agg(lambda values: ";".join(sorted(map(str, values))))
        .rename("passed_biological_sample_ids")
        .reset_index()
    )
    recurrence = recurrence.merge(
        passed_ids, on=RECURRENCE_KEY, how="left", validate="one_to_one"
    )
    recurrence["passed_biological_sample_ids"] = recurrence[
        "passed_biological_sample_ids"
    ].fillna("")
    recurrence["aggregation"] = np.where(
        recurrence["resource_axis"].eq(1),
        "all_5_technical_replicates",
        "not_applicable_fixed_fake",
    )
    score = score.merge(
        recurrence[
            RECURRENCE_KEY
            + [
                "biological_samples_total",
                "biological_samples_joint_pass_all_replicates",
                "joint_pair_null_recurrence_fraction",
                "passed_biological_sample_ids",
                "aggregation",
            ]
        ],
        on=RECURRENCE_KEY,
        how="left",
        validate="many_to_one",
    ).rename(
        columns={
            "biological_samples_joint_pass_all_replicates": (
                "joint_pair_null_pass_biological_samples"
            )
        }
    )
    score["score_eligible"] = score["resource_axis"].astype(np.int8)
    score[SCORE] = (
        0.50 * score["pair_null_support_fraction"]
        + 0.20 * score["method_support_fraction_all"]
        + 0.20 * score["joint_pair_null_recurrence_fraction"]
        + 0.10 * score["global_null_support_fraction"]
    )
    high = (
        score["score_eligible"].eq(1)
        & score[SCORE].ge(0.75)
        & score["observed_local_lr_score"].gt(0)
        & score["pair_null_pass_count"].eq(3)
        & score["methods_detected"].ge(2)
        & score["joint_pair_null_pass_biological_samples"].ge(2)
    )
    medium = (
        score["score_eligible"].eq(1)
        & ~high
        & score[SCORE].ge(0.50)
        & score["observed_local_lr_score"].gt(0)
        & score["pair_null_pass_count"].ge(2)
        & score["methods_detected"].ge(1)
    )
    score[TIER] = np.select(
        [score["score_eligible"].eq(0), high, medium],
        [
            "out_of_scope_fixed_fake",
            "high_f5_full",
            "medium_f5_full",
        ],
        default="low_f5_full",
    )
    score.loc[score["score_eligible"].eq(0), SCORE] = np.nan
    score = score.sort_values(ROW_KEY).reset_index(drop=True)

    metric, sample_metric, bootstrap = build_metrics(score)
    topk = build_topk(score)
    fdp, hard_negative, condition_performance = build_selection_audits(score)
    monotonicity, monotonicity_summary = build_monotonicity(score)
    rank_detail, rank_summary = build_rank_stability(score)
    fixed = score[score["resource_axis"].eq(0)]
    fixed_summary = pd.DataFrame(
        [
            {
                "rows": len(fixed),
                "samples": fixed["sample_id"].nunique(),
                "technical_replicates": fixed[
                    ["sample_id", "replicate"]
                ].drop_duplicates().shape[0],
                "method_detections": int(fixed["methods_detected"].sum()),
                "score_values_present": int(fixed[SCORE].notna().sum()),
                "tier": "out_of_scope_fixed_fake",
                "primary_metrics_excluded": 1,
            }
        ]
    )

    output.mkdir(parents=True, exist_ok=True)
    recurrence.to_csv(output / "biological_recurrence.tsv", sep="\t", index=False)
    technical.to_csv(
        output / "technical_replicate_recurrence.tsv", sep="\t", index=False
    )
    metric.to_csv(output / "metric_summary.tsv", sep="\t", index=False)
    sample_metric.to_csv(output / "metric_by_sample.tsv", sep="\t", index=False)
    bootstrap.to_csv(output / "section_bootstrap_metrics.tsv", sep="\t", index=False)
    topk.to_csv(output / "topk_by_sample_replicate.tsv", sep="\t", index=False)
    fdp.to_csv(output / "empirical_fdp.tsv", sep="\t", index=False)
    hard_negative.to_csv(
        output / "hard_negative_false_positive_rates.tsv", sep="\t", index=False
    )
    condition_performance.to_csv(
        output / "condition_selection_rates.tsv", sep="\t", index=False
    )
    monotonicity.to_csv(
        output / "perturbation_monotonicity.tsv", sep="\t", index=False
    )
    monotonicity_summary.to_csv(
        output / "perturbation_monotonicity_summary.tsv", sep="\t", index=False
    )
    rank_detail.to_csv(
        output / "technical_replicate_rank_stability.tsv", sep="\t", index=False
    )
    rank_summary.to_csv(
        output / "technical_replicate_rank_stability_summary.tsv",
        sep="\t",
        index=False,
    )
    fixed_summary.to_csv(
        output / "fixed_fake_out_of_scope_summary.tsv", sep="\t", index=False
    )

    output_score = Path(args.output_score).resolve()
    output_meta = Path(args.output_meta).resolve()
    output_score.parent.mkdir(parents=True, exist_ok=True)
    score.to_csv(output_score, sep="\t", index=False)
    tier_counts = score[TIER].value_counts().to_dict()
    summary = {
        "status": "PASS_F5_SEMISIM_FULL_SCORE_BUILT",
        "run_id": args.run_id,
        "truth_axis_rows": len(score),
        "resource_score_rows": int(score["score_eligible"].sum()),
        "fixed_fake_out_of_scope_rows": int(score["score_eligible"].eq(0).sum()),
        "biological_samples": score["sample_id"].nunique(),
        "technical_replicates_per_sample": 5,
        "conditions": score[
            ["sample_id", "replicate", "condition_id"]
        ].drop_duplicates().shape[0],
        "recurrence_keys": len(recurrence),
        "tier_counts": {key: int(value) for key, value in tier_counts.items()},
        "metric_rows": len(metric),
        "bootstrap_rows": len(bootstrap),
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "validation_pending": True,
        "boundary": (
            "Frozen computational score on truth-labelled semi-simulation; "
            "not a calibrated probability or experimental validation."
        ),
    }
    (output / "build_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    metadata = {
        **summary,
        "status": "f5_semisim_full_score",
        "formula": (
            "0.50*pair_null_support_fraction + "
            "0.20*method_support_fraction_all + "
            "0.20*joint_pair_null_recurrence_fraction + "
            "0.10*global_null_support_fraction"
        ),
        "recurrence_denominator": 3,
        "recurrence_aggregation": "all_5_technical_replicates",
        "pair_fdr_threshold": 0.05,
        "global_fdr_threshold": 0.05,
        "top_k": TOP_K,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_cluster": "sample_id",
        "fixed_fake_policy": (
            "Retained as out-of-scope negative control; score is NA and rows "
            "are excluded from primary AUROC/AUPRC."
        ),
        "probability_mapping": None,
        "brier_score_computed": False,
        "source_sha256": {
            "method_audit": digest(method_path),
            "method_validation_summary": digest(
                Path(args.method_validation_summary)
            ),
            "null_audit": digest(null_audit_path),
            "null_validation_summary": digest(
                Path(args.null_validation_summary)
            ),
            "tf6_development_score": digest(
                Path(args.tf6_development_score)
            ),
            "v3_baseline_checksums": digest(Path(args.v3_checksums)),
        },
    }
    output_meta.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
