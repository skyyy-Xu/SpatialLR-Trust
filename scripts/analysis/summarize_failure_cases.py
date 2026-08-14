#!/usr/bin/env python3
"""Summarize SpatialLR-Trust filtered and failure-case diagnostics."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "20260709_2230_failure-cases")
OUT_DIR = PROJECT / "results/task_e_failure_cases"
FIG_DIR = PROJECT / "figures/task_e_failure_cases"
RUN_DIR = PROJECT / "runs" / RUN_ID
BASE = "spatiallr_failure_cases"

PALETTE = {
    "blue": "#4C78A8",
    "green": "#59A14F",
    "red": "#E15759",
    "orange": "#F28E2B",
    "purple": "#B07AA1",
    "gray": "#B9B9B9",
    "dark": "#444444",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(PROJECT / rel, sep="\t")


def clean_label(text: str) -> str:
    return (
        text.replace("both_external_supported_but_failed_null_gate", "External support, null failure")
        .replace("failed_all_three_null_gates", "Failed all null gates")
        .replace("raw_score_supported_but_no_external_method_support", "Raw score, no external support")
        .replace("v2_high_downgraded_without_both_external_methods", "v2 high downgraded")
        .replace("low_percentile_and_low_enrichment", "Low percentile + enrichment")
        .replace("low_receiver_enrichment", "Low receiver enrichment")
        .replace("low_receiver_percentile", "Low receiver percentile")
        .replace("target_supported", "Target-supported high v3")
    )


def candidate_name(row: pd.Series) -> str:
    return f"{row['sender_compartment']}->{row['receiver_compartment']} {row['ligand']}-{row['receptor']}"


def build_outputs() -> dict[str, object]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    filtered = read_tsv("results/task_e_score_diagnostics/v3_filtered_case_examples.tsv")
    filter_summary = read_tsv("results/task_e_score_diagnostics/v3_filtered_case_type_summary.tsv")
    gaps = read_tsv("results/task_e_target_activation/target_activation_high_v3_gap_candidates.tsv")
    gap_summary = read_tsv("results/task_e_target_activation/target_activation_high_v3_gap_reason_summary.tsv")
    scored = read_tsv("results/task_e_target_activation/target_activation_all_score_v3_with_target.tsv")

    filtered["score_v3"] = pd.to_numeric(filtered["score_v3"], errors="coerce")
    filtered["null_pass_count"] = pd.to_numeric(filtered["null_pass_count"], errors="coerce")
    gaps["spatiallr_trust_score_pilot_v3"] = pd.to_numeric(gaps["spatiallr_trust_score_pilot_v3"], errors="coerce")
    gaps["receiver_target_enrichment"] = pd.to_numeric(gaps["receiver_target_enrichment"], errors="coerce")
    gaps["receiver_target_percentile"] = pd.to_numeric(gaps["receiver_target_percentile"], errors="coerce")

    filter_summary = filter_summary.copy()
    filter_summary["filter_case_label"] = filter_summary["filter_case_type"].map(clean_label)
    filter_summary["fraction_of_scored_candidates"] = filter_summary["candidate_rows"] / len(scored)
    filter_summary["fraction_of_filtered_examples"] = filter_summary["candidate_rows"] / len(filtered)
    filter_summary.to_csv(OUT_DIR / "failure_case_type_summary.tsv", sep="\t", index=False)

    top_filtered = (
        filtered.sort_values(["filter_case_type", "score_v3"], ascending=[True, False])
        .groupby("filter_case_type", as_index=False)
        .head(5)
        .copy()
    )
    top_filtered["candidate"] = top_filtered.apply(candidate_name, axis=1)
    top_filtered["filter_case_label"] = top_filtered["filter_case_type"].map(clean_label)
    top_filtered[
        [
            "filter_case_label", "filter_case_type", "candidate", "dataset", "sample_id", "cancer",
            "pathway", "score_v3", "tier_v3", "null_pass_count", "null_failure_reason",
            "commot_supported", "liana_supported", "external_method_support_count",
            "dataset_recurrence_n", "evidence_label",
        ]
    ].to_csv(OUT_DIR / "failure_case_top_examples.tsv", sep="\t", index=False)

    gap_summary = gap_summary.copy()
    gap_summary["target_gap_label"] = gap_summary["target_gap_reason"].map(clean_label)
    gap_summary["fraction_of_high_v3"] = gap_summary["rows"] / gap_summary["rows"].sum()
    gap_summary.to_csv(OUT_DIR / "high_v3_target_gap_summary.tsv", sep="\t", index=False)

    top_gaps = (
        gaps.sort_values(["target_gap_reason", "spatiallr_trust_score_pilot_v3"], ascending=[True, False])
        .groupby("target_gap_reason", as_index=False)
        .head(5)
        .copy()
    )
    top_gaps["target_gap_label"] = top_gaps["target_gap_reason"].map(clean_label)
    top_gaps[
        [
            "target_gap_label", "target_gap_reason", "candidate", "dataset", "sample_id", "cancer",
            "pathway", "spatiallr_trust_score_pilot_v3", "tri_method_support_count",
            "external_method_support_count", "spatial_null_p", "label_null_p", "fake_lr_null_p",
            "dataset_recurrence_n", "target_genes_detected", "receiver_target_enrichment",
            "receiver_target_percentile", "target_activation_support",
        ]
    ].to_csv(OUT_DIR / "high_v3_target_gap_examples.tsv", sep="\t", index=False)

    null_dist = (
        filtered.groupby(["filter_case_type", "null_pass_count"], dropna=False)
        .size()
        .reset_index(name="candidate_rows")
    )
    null_dist["filter_case_label"] = null_dist["filter_case_type"].map(clean_label)
    null_dist.to_csv(OUT_DIR / "failure_case_null_pass_distribution.tsv", sep="\t", index=False)

    high_count = int((scored["confidence_tier_pilot_v3"] == "high_pilot_v3").sum())
    target_supported_high = int(gap_summary.loc[gap_summary["target_gap_reason"] == "target_supported", "rows"].sum())
    unsupported_high = high_count - target_supported_high
    metrics = pd.DataFrame([
        {"metric": "scored_candidate_rows", "value": len(scored)},
        {"metric": "filtered_case_examples", "value": len(filtered)},
        {"metric": "filter_case_types", "value": filter_summary.shape[0]},
        {"metric": "external_supported_but_failed_null_gate", "value": int(filter_summary.loc[filter_summary["filter_case_type"] == "both_external_supported_but_failed_null_gate", "candidate_rows"].sum())},
        {"metric": "failed_all_three_null_gates", "value": int(filter_summary.loc[filter_summary["filter_case_type"] == "failed_all_three_null_gates", "candidate_rows"].sum())},
        {"metric": "raw_score_supported_but_no_external_method_support", "value": int(filter_summary.loc[filter_summary["filter_case_type"] == "raw_score_supported_but_no_external_method_support", "candidate_rows"].sum())},
        {"metric": "v2_high_downgraded_without_both_external_methods", "value": int(filter_summary.loc[filter_summary["filter_case_type"] == "v2_high_downgraded_without_both_external_methods", "candidate_rows"].sum())},
        {"metric": "high_v3_candidate_rows", "value": high_count},
        {"metric": "high_v3_target_supported", "value": target_supported_high},
        {"metric": "high_v3_target_unsupported", "value": unsupported_high},
    ])
    metrics.to_csv(OUT_DIR / "failure_case_key_metrics.tsv", sep="\t", index=False)

    draw_figure(filter_summary, null_dist, gap_summary, top_filtered, top_gaps)

    manifest = {
        "run_id": RUN_ID,
        "outputs": {
            "tables": [str(p.relative_to(PROJECT)) for p in sorted(OUT_DIR.glob("*.tsv"))],
            "figures": [str(p.relative_to(PROJECT)) for p in sorted(FIG_DIR.glob(f"{BASE}.*"))],
        },
        "key_metrics": {row["metric"]: int(row["value"]) for _, row in metrics.iterrows()},
        "interpretation_boundary": "Failure cases explain computational filtering and target-proxy gaps; they do not prove absence of biological communication.",
    }
    (FIG_DIR / f"{BASE}_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def style_barh(ax) -> None:
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=6, length=2)


def add_label(ax, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=10, fontweight="bold")


def draw_figure(filter_summary: pd.DataFrame, null_dist: pd.DataFrame, gap_summary: pd.DataFrame, top_filtered: pd.DataFrame, top_gaps: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(8.1, 6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    fs = filter_summary.sort_values("candidate_rows", ascending=True)
    colors = [PALETTE["purple"], PALETTE["orange"], PALETTE["red"], PALETTE["blue"]][-len(fs):]
    bars = ax_a.barh(fs["filter_case_label"], fs["candidate_rows"], color=colors, height=0.62)
    ax_a.set_xlabel("Candidate rows", fontsize=6.5)
    ax_a.set_title("Filtered candidates expose distinct failure modes", loc="left", fontsize=8)
    ax_a.set_xscale("log")
    style_barh(ax_a)
    for bar, val in zip(bars, fs["candidate_rows"]):
        ax_a.text(float(val) * 1.05, bar.get_y() + bar.get_height() / 2, f"{int(val):,}", va="center", fontsize=6)
    add_label(ax_a, "a")

    pivot = null_dist.pivot_table(index="filter_case_label", columns="null_pass_count", values="candidate_rows", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(fs["filter_case_label"])
    bottom = pd.Series(0, index=pivot.index, dtype=float)
    null_colors = [PALETTE["red"], PALETTE["orange"], PALETTE["gray"], PALETTE["green"]]
    for i, col in enumerate(sorted(pivot.columns)):
        vals = pivot[col]
        ax_b.barh(pivot.index, vals, left=bottom, color=null_colors[int(col) if int(col) < len(null_colors) else -1], height=0.62, label=f"{int(col)} gates")
        bottom += vals
    ax_b.set_xlabel("Candidate rows", fontsize=6.5)
    ax_b.set_title("Null gates separate support from instability", loc="left", fontsize=8)
    ax_b.set_xscale("log")
    style_barh(ax_b)
    ax_b.legend(title="Null gates passed", title_fontsize=5.8, fontsize=5.8, loc="lower right")
    add_label(ax_b, "b")

    gs2 = gap_summary.copy()
    gs2 = gs2.sort_values("rows", ascending=True)
    gap_colors = [PALETTE["orange"] if r != "target_supported" else PALETTE["green"] for r in gs2["target_gap_reason"]]
    bars = ax_c.barh(gs2["target_gap_label"], gs2["rows"], color=gap_colors, height=0.62)
    ax_c.set_xlabel("High-v3 candidate rows", fontsize=6.5)
    ax_c.set_title("Target proxy reveals gaps among high-scoring candidates", loc="left", fontsize=8)
    style_barh(ax_c)
    for bar, val in zip(bars, gs2["rows"]):
        ax_c.text(float(val) + 2, bar.get_y() + bar.get_height() / 2, f"{int(val):,}", va="center", fontsize=6)
    add_label(ax_c, "c")

    ax_d.axis("off")
    ax_d.set_title("Representative filtered candidates remain inspectable", loc="left", fontsize=8, pad=2)
    reps = []
    selected_filter_types = [
        "both_external_supported_but_failed_null_gate",
        "failed_all_three_null_gates",
        "raw_score_supported_but_no_external_method_support",
    ]
    for ftype in selected_filter_types:
        subset = top_filtered[top_filtered["filter_case_type"] == ftype].sort_values("score_v3", ascending=False)
        if not subset.empty:
            row = subset.iloc[0]
            if ftype == "failed_all_three_null_gates":
                reason = "spatial, label and fake-LR gates failed"
            elif ftype == "raw_score_supported_but_no_external_method_support":
                reason = "no external method support"
            else:
                reason = str(row["null_failure_reason"]).replace("fake_lr_null_p", "fake-LR null p")
            reps.append([
                clean_label(ftype),
                row["candidate"].replace("stroma_fibroblast", "stroma"),
                reason,
                f"{row['score_v3']:.3f}",
            ])
    gap_subset = top_gaps[top_gaps["target_gap_reason"] != "target_supported"].sort_values("spatiallr_trust_score_pilot_v3", ascending=False)
    if not gap_subset.empty:
        row = gap_subset.iloc[0]
        reps.append([
            clean_label(row["target_gap_reason"]),
            row["candidate"].replace("stroma_fibroblast", "stroma"),
            f"target enrich {row['receiver_target_enrichment']:.2f}",
            f"{row['spatiallr_trust_score_pilot_v3']:.3f}",
        ])
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    row_colors = [PALETTE["orange"], PALETTE["red"], PALETTE["blue"], PALETTE["purple"]]
    for index, (mode, candidate, reason, score) in enumerate(reps):
        y = 0.88 - index * 0.21
        ax_d.text(0.02, y, mode, fontsize=5.9, fontweight="bold", color=row_colors[index], va="top")
        ax_d.text(0.02, y - 0.065, candidate, fontsize=5.5, color=PALETTE["dark"], va="top")
        ax_d.text(0.02, y - 0.125, f"{reason}; v3 = {score}", fontsize=5.3, color="#666666", va="top")
        if index < len(reps) - 1:
            ax_d.plot([0.02, 0.98], [y - 0.175, y - 0.175], color="#E3E3E3", linewidth=0.5)
    ax_d.text(0.02, 0.01, "Examples explain filtering; they do not prove biological absence.", transform=ax_d.transAxes, fontsize=5.4, color="#555555")
    add_label(ax_d, "d")

    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIG_DIR / f"{BASE}.{ext}", dpi=320, bbox_inches="tight")
    try:
        fig.savefig(FIG_DIR / f"{BASE}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(FIG_DIR / f"{BASE}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    manifest = build_outputs()
    print(json.dumps(manifest, indent=2))
