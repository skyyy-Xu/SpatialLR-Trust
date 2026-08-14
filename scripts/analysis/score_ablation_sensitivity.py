#!/usr/bin/env python3
"""Score-component ablation analysis for SpatialLR-Trust pilot v3."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "20260709_2305_score-ablation")
INPUT = PROJECT / "results/task_e_target_activation/target_activation_all_score_v3_with_target.tsv"
OUT_DIR = PROJECT / "results/task_e_score_ablation"
FIG_DIR = PROJECT / "figures/task_e_score_ablation"
RUN_DIR = PROJECT / "runs" / RUN_ID
BASE = "spatiallr_score_ablation"

WEIGHTS = {
    "null_support_mean": 0.40,
    "pilot_score_percentile": 0.15,
    "recurrence_support": 0.15,
    "annotation_quality": 0.10,
    "external_method_support_fraction": 0.20,
}
LABELS = {
    "full_v3": "Full v3",
    "no_null": "No null layer",
    "no_external": "No external support",
    "no_recurrence": "No recurrence",
    "no_raw_score": "No raw score",
    "no_annotation": "No annotation quality",
    "equal_weight": "Equal weights",
}
PALETTE = {
    "blue": "#4C78A8",
    "light_blue": "#9ECAE1",
    "green": "#59A14F",
    "orange": "#F28E2B",
    "red": "#E15759",
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
    "legend.frameon": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def candidate_id(df: pd.DataFrame) -> pd.Series:
    fields = ["dataset", "sample_id", "sender_compartment", "receiver_compartment", "ligand", "receptor", "pathway"]
    return df[fields].astype(str).agg("|".join, axis=1)


def renorm_score(df: pd.DataFrame, excluded: set[str]) -> pd.Series:
    cols = [c for c in WEIGHTS if c not in excluded]
    denom = sum(WEIGHTS[c] for c in cols)
    return sum(df[c] * WEIGHTS[c] for c in cols) / denom


def build_tables() -> tuple[pd.DataFrame, dict[str, object]]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t")
    for col in list(WEIGHTS) + [
        "spatiallr_trust_score_pilot_v3", "spatial_null_p", "label_null_p", "fake_lr_null_p",
        "target_activation_support", "external_method_support_count", "receiver_target_enrichment",
        "receiver_target_percentile",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df["candidate_id"] = candidate_id(df)
    df["all_null_gates"] = (df["spatial_null_p"] <= 0.10) & (df["label_null_p"] <= 0.10) & (df["fake_lr_null_p"] <= 0.10)
    df["both_external"] = df["external_method_support_count"] == 2
    df["full_high"] = df["confidence_tier_pilot_v3"] == "high_pilot_v3"

    df["full_v3"] = df["spatiallr_trust_score_pilot_v3"]
    df["no_null"] = renorm_score(df, {"null_support_mean"})
    df["no_external"] = renorm_score(df, {"external_method_support_fraction"})
    df["no_recurrence"] = renorm_score(df, {"recurrence_support"})
    df["no_raw_score"] = renorm_score(df, {"pilot_score_percentile"})
    df["no_annotation"] = renorm_score(df, {"annotation_quality"})
    df["equal_weight"] = df[list(WEIGHTS)].mean(axis=1)

    full_high_ids = set(df.loc[df["full_high"], "candidate_id"])
    n_high = len(full_high_ids)
    variants = ["full_v3", "no_null", "no_external", "no_recurrence", "no_raw_score", "no_annotation", "equal_weight"]
    rows = []
    top_rows = []
    for var in variants:
        top = df.sort_values(var, ascending=False).head(n_high).copy()
        top_ids = set(top["candidate_id"])
        overlap = len(full_high_ids & top_ids)
        rows.append({
            "variant": var,
            "variant_label": LABELS[var],
            "top_n": n_high,
            "overlap_with_full_high": overlap,
            "retention_fraction": overlap / n_high if n_high else 0,
            "spearman_with_full_v3": df["full_v3"].corr(df[var], method="spearman"),
            "top_n_all_null_gate_rows": int(top["all_null_gates"].sum()),
            "top_n_both_external_rows": int(top["both_external"].sum()),
            "top_n_target_supported_rows": int(top["target_activation_support"].sum()),
            "top_n_target_supported_fraction": float(top["target_activation_support"].mean()),
            "top_n_mean_full_v3_score": float(top["full_v3"].mean()),
            "top_n_mean_variant_score": float(top[var].mean()),
        })
        keep_cols = [
            "candidate_id", "dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment",
            "ligand", "receptor", "pathway", "full_v3", var, "confidence_tier_pilot_v3", "all_null_gates",
            "both_external", "target_activation_support",
        ]
        top["variant"] = var
        top["rank_in_variant"] = np.arange(1, len(top) + 1)
        ordered_cols = []
        for col in ["variant", "rank_in_variant"] + keep_cols:
            if col not in ordered_cols:
                ordered_cols.append(col)
        top_rows.append(top[ordered_cols])
    ablation = pd.DataFrame(rows)
    ablation.to_csv(OUT_DIR / "score_ablation_topn_overlap.tsv", sep="\t", index=False)
    pd.concat(top_rows, ignore_index=True).to_csv(OUT_DIR / "score_ablation_topn_candidates.tsv", sep="\t", index=False)

    gate_rows = []
    gate_defs = [
        ("full_high_rule", "Full high rule", (df["full_v3"] >= 0.75) & df["all_null_gates"] & df["both_external"]),
        ("score_threshold_only", "Score >= 0.75 only", df["full_v3"] >= 0.75),
        ("drop_null_gate", "Drop null gate", (df["full_v3"] >= 0.75) & df["both_external"]),
        ("drop_external_gate", "Drop external gate", (df["full_v3"] >= 0.75) & df["all_null_gates"]),
        ("drop_both_gates", "Drop both gates", df["full_v3"] >= 0.75),
    ]
    for key, label, mask in gate_defs:
        subset = df[mask].copy()
        gate_rows.append({
            "gate_scenario": key,
            "gate_label": label,
            "candidate_rows": len(subset),
            "target_supported_rows": int(subset["target_activation_support"].sum()) if len(subset) else 0,
            "target_supported_fraction": float(subset["target_activation_support"].mean()) if len(subset) else 0,
            "all_null_gate_rows": int(subset["all_null_gates"].sum()) if len(subset) else 0,
            "both_external_rows": int(subset["both_external"].sum()) if len(subset) else 0,
        })
    gates = pd.DataFrame(gate_rows)
    gates.to_csv(OUT_DIR / "score_gate_sensitivity.tsv", sep="\t", index=False)

    component_rows = []
    tiers = ["high_pilot_v3", "medium_pilot_v3", "low_pilot_v3"]
    for tier in tiers:
        subset = df[df["confidence_tier_pilot_v3"] == tier]
        for comp, weight in WEIGHTS.items():
            component_rows.append({
                "tier": tier,
                "component": comp,
                "component_label": comp.replace("_", " "),
                "weight": weight,
                "mean_value": float(subset[comp].mean()),
                "median_value": float(subset[comp].median()),
                "rows": len(subset),
            })
    components = pd.DataFrame(component_rows)
    components.to_csv(OUT_DIR / "score_component_by_tier.tsv", sep="\t", index=False)

    rank_df = df[["candidate_id", "dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment", "ligand", "receptor", "pathway", "full_high"] + variants].copy()
    for var in variants:
        rank_df[f"rank_{var}"] = rank_df[var].rank(method="first", ascending=False).astype(int)
    rank_df["rank_shift_no_null_vs_full"] = rank_df["rank_no_null"] - rank_df["rank_full_v3"]
    rank_df["rank_shift_no_external_vs_full"] = rank_df["rank_no_external"] - rank_df["rank_full_v3"]
    rank_df.to_csv(OUT_DIR / "score_ablation_candidate_ranks.tsv", sep="\t", index=False)

    metrics = {
        "run_id": RUN_ID,
        "input": str(INPUT.relative_to(PROJECT)),
        "candidate_rows": int(len(df)),
        "full_high_rows": int(n_high),
        "full_high_target_supported_rows": int(df.loc[df["full_high"], "target_activation_support"].sum()),
        "score_threshold_only_rows": int(gates.loc[gates["gate_scenario"] == "score_threshold_only", "candidate_rows"].iloc[0]),
        "drop_null_gate_rows": int(gates.loc[gates["gate_scenario"] == "drop_null_gate", "candidate_rows"].iloc[0]),
        "drop_external_gate_rows": int(gates.loc[gates["gate_scenario"] == "drop_external_gate", "candidate_rows"].iloc[0]),
        "no_null_topn_retention": float(ablation.loc[ablation["variant"] == "no_null", "retention_fraction"].iloc[0]),
        "no_external_topn_retention": float(ablation.loc[ablation["variant"] == "no_external", "retention_fraction"].iloc[0]),
        "lowest_retention_variant": str(ablation.sort_values("retention_fraction").iloc[0]["variant"]),
        "lowest_retention_fraction": float(ablation.sort_values("retention_fraction").iloc[0]["retention_fraction"]),
        "outputs": {
            "ablation_summary": str((OUT_DIR / "score_ablation_topn_overlap.tsv").relative_to(PROJECT)),
            "gate_sensitivity": str((OUT_DIR / "score_gate_sensitivity.tsv").relative_to(PROJECT)),
            "component_by_tier": str((OUT_DIR / "score_component_by_tier.tsv").relative_to(PROJECT)),
            "candidate_ranks": str((OUT_DIR / "score_ablation_candidate_ranks.tsv").relative_to(PROJECT)),
        },
        "boundary": "Ablations test score-component sensitivity, not calibrated causal contributions or ground-truth accuracy.",
    }
    (OUT_DIR / "score_ablation_summary.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    return df, ablation, gates, components, metrics


def add_label(ax, label: str) -> None:
    ax.text(-0.11, 1.05, label, transform=ax.transAxes, fontsize=10, fontweight="bold", ha="left", va="bottom")


def style(ax) -> None:
    ax.tick_params(axis="both", labelsize=6, length=2)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    ax.set_axisbelow(True)


def draw_figure(ablation: pd.DataFrame, gates: pd.DataFrame, components: pd.DataFrame, metrics: dict[str, object]) -> None:
    fig = plt.figure(figsize=(8.3, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_ab = ablation[ablation["variant"] != "full_v3"].sort_values("retention_fraction")
    bars = ax_a.barh(plot_ab["variant_label"], plot_ab["retention_fraction"], color=PALETTE["blue"], height=0.62)
    ax_a.set_xlim(0, 1.05)
    ax_a.set_xlabel("Top-N retention vs full high-v3 tier", fontsize=6.5)
    ax_a.set_title("Removing score layers changes the high-priority set", loc="left", fontsize=8)
    style(ax_a)
    for bar, val in zip(bars, plot_ab["retention_fraction"]):
        ax_a.text(val + 0.02, bar.get_y() + bar.get_height() / 2, f"{val:.2f}", va="center", fontsize=6)
    add_label(ax_a, "a")

    gate_order = ["full_high_rule", "drop_external_gate", "drop_null_gate", "score_threshold_only"]
    g = gates.set_index("gate_scenario").loc[gate_order].reset_index()
    unsupported = g["candidate_rows"] - g["target_supported_rows"]
    y = np.arange(len(g))
    ax_b.barh(y, g["target_supported_rows"], color=PALETTE["green"], label="Target-supported", height=0.62)
    ax_b.barh(y, unsupported, left=g["target_supported_rows"], color=PALETTE["gray"], label="No target proxy", height=0.62)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(g["gate_label"], fontsize=6)
    ax_b.invert_yaxis()
    ax_b.set_xlabel("Candidate rows", fontsize=6.5)
    ax_b.set_title("Gate relaxation expands candidate counts", loc="left", fontsize=8)
    style(ax_b)
    ax_b.legend(fontsize=5.8, loc="lower right")
    for i, val in enumerate(g["candidate_rows"]):
        ax_b.text(val + 8, i, f"{int(val):,}", va="center", fontsize=6)
    add_label(ax_b, "b")

    comp_order = ["null_support_mean", "external_method_support_fraction", "recurrence_support", "pilot_score_percentile"]
    comp_labels = ["Null", "External", "Recurrence", "Raw percentile"]
    pivot = components[components["component"].isin(comp_order)].pivot(index="component", columns="tier", values="mean_value").loc[comp_order]
    x = np.arange(len(comp_order))
    width = 0.24
    colors = [PALETTE["green"], PALETTE["orange"], PALETTE["gray"]]
    for i, tier in enumerate(["high_pilot_v3", "medium_pilot_v3", "low_pilot_v3"]):
        ax_c.bar(x + (i - 1) * width, pivot[tier], width=width, color=colors[i], label=tier.replace("_pilot_v3", ""))
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(comp_labels, rotation=25, ha="right", fontsize=6)
    ax_c.set_ylim(0, 1.05)
    ax_c.set_ylabel("Mean component value", fontsize=6.5)
    ax_c.set_title("High-v3 candidates are enriched across evidence layers", loc="left", fontsize=8)
    ax_c.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    ax_c.legend(fontsize=5.8, loc="upper right")
    add_label(ax_c, "c")

    ax_d.axis("off")
    ax_d.set_title("Sensitivity summary", loc="left", fontsize=8)
    rows = [
        ["Full high-v3 rows", f"{metrics['full_high_rows']:,}"],
        ["Target-supported high-v3", f"{metrics['full_high_target_supported_rows']:,}"],
        ["Score >= 0.75 only", f"{metrics['score_threshold_only_rows']:,}"],
        ["Drop null gate", f"{metrics['drop_null_gate_rows']:,}"],
        ["Drop external gate", f"{metrics['drop_external_gate_rows']:,}"],
        ["Lowest top-N retention", f"{metrics['lowest_retention_variant']} ({metrics['lowest_retention_fraction']:.2f})"],
    ]
    table = ax_d.table(cellText=rows, colLabels=["Metric", "Value"], cellLoc="left", colLoc="left", loc="upper left", bbox=[0, 0.18, 0.98, 0.72], colWidths=[0.62, 0.36])
    table.auto_set_font_size(False)
    table.set_fontsize(6.2)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        cell.set_linewidth(0.4)
        if row == 0:
            cell.set_facecolor("#F2F2F2")
            cell.set_text_props(weight="bold")
    ax_d.text(0, 0.05, "Ablations test score sensitivity, not calibrated causal contributions or ground-truth accuracy.", transform=ax_d.transAxes, fontsize=5.8, color=PALETTE["dark"])
    add_label(ax_d, "d")

    for ext in ["svg", "pdf", "png"]:
        fig.savefig(FIG_DIR / f"{BASE}.{ext}", dpi=320, bbox_inches="tight")
    try:
        fig.savefig(FIG_DIR / f"{BASE}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    except TypeError:
        fig.savefig(FIG_DIR / f"{BASE}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df, ablation, gates, components, metrics = build_tables()
    draw_figure(ablation, gates, components, metrics)
    manifest = {
        "run_id": RUN_ID,
        "figure_contract": {
            "core_conclusion": "SpatialLR-Trust high-priority selection depends most on null and external-support layers; gate relaxation expands candidate counts.",
            "archetype": "quantitative grid",
            "backend": "python/matplotlib",
            "boundary": metrics["boundary"],
        },
        "outputs": {
            "tables": [str(p.relative_to(PROJECT)) for p in sorted(OUT_DIR.glob("*.tsv"))] + [str((OUT_DIR / "score_ablation_summary.json").relative_to(PROJECT))],
            "figures": [str(p.relative_to(PROJECT)) for p in sorted(FIG_DIR.glob(f"{BASE}.*"))],
        },
        "key_metrics": {k: v for k, v in metrics.items() if not isinstance(v, dict)},
    }
    (FIG_DIR / f"{BASE}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
