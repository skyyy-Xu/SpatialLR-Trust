#!/usr/bin/env python3
"""Render the global SpatialLR-Trust Results overview figure from derived outputs."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "20260709_2340_global-results-figure")
OUT_DIR = PROJECT / "results/task_e_global_results_overview"
FIG_DIR = PROJECT / "figures/task_e_global_results_overview"
RUN_DIR = PROJECT / "runs" / RUN_ID
BASE_NAME = "spatiallr_global_results_overview"

PALETTE = {
    "blue": "#4C78A8",
    "light_blue": "#9ECAE1",
    "green": "#59A14F",
    "light_green": "#A1D99B",
    "red": "#E15759",
    "orange": "#F28E2B",
    "yellow": "#EDC948",
    "purple": "#B07AA1",
    "gray": "#B9B9B9",
    "dark_gray": "#555555",
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


def read_json(rel: str) -> dict:
    with (PROJECT / rel).open() as handle:
        return json.load(handle)


def read_tsv(rel: str) -> list[dict[str, str]]:
    with (PROJECT / rel).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_metrics() -> dict:
    pilot = read_json("results/task_c_pilot_baseline/stdlib_marker_lr_pilot_summary.json")
    commot = read_json("results/task_c_commot_baseline/commot_all_summary.json")
    liana = read_json("results/task_c_liana_baseline/liana_cellphonedb_all_summary.json")
    tri = read_json("results/task_c_method_consistency/tri_method_consistency_summary.json")
    spatial_null = read_json("runs/20260709_0115_spatial-null-pilot/spatial_null_summary.json")
    label_null = read_json("runs/20260709_1137_label-null-pilot/label_null_summary.json")
    fake_null = read_json("runs/20260709_1511_fake-lr-null-pilot/fake_lr_null_summary.json")
    score = read_json("results/task_e_scores/spatiallr_trust_score_pilot_v3_summary.json")
    diag = read_json("results/task_e_score_diagnostics/v3_benchmark_diagnostics_summary.json")
    target = read_json("results/task_e_target_activation/target_activation_all_summary.json")
    gap = read_json("results/task_e_target_activation/target_activation_high_v3_gap_summary.json")
    public = read_json("results/task_e_target_activation/public_evidence_top40_summary.json")
    case = read_json("results/task_e_target_activation/same_cancer_case_study_summary.json")
    audit = read_json("results/task_e_target_activation/same_cancer_original_paper_evidence_summary.json")
    qc_rows = read_tsv("docs/minimal_input_qc.tsv")
    return {
        "pilot": pilot,
        "commot": commot,
        "liana": liana,
        "tri": tri,
        "spatial_null": spatial_null,
        "label_null": label_null,
        "fake_null": fake_null,
        "score": score,
        "diag": diag,
        "target": target,
        "gap": gap,
        "public": public,
        "case": case,
        "audit": audit,
        "qc_rows": qc_rows,
    }


def build_source_tables(m: dict) -> dict[str, Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tables: dict[str, Path] = {}

    dataset_counts: dict[str, int] = {}
    for row in m["qc_rows"]:
        dataset_counts[row["dataset"]] = dataset_counts.get(row["dataset"], 0) + 1
    rows = [{"dataset": k, "samples": v} for k, v in sorted(dataset_counts.items())]
    path = OUT_DIR / "global_overview_dataset_samples.tsv"
    write_tsv(path, rows, ["dataset", "samples"])
    tables["dataset_samples"] = path

    rows = [
        {"stage": "Marker baseline", "count": m["pilot"]["candidates"], "note": "scored candidate universe"},
        {"stage": "All null gates", "count": m["diag"]["null_all_pass_rows"], "note": "passed current null gates"},
        {"stage": "High v3", "count": m["score"]["tiers"]["high_pilot_v3"], "note": "operational high-priority v3 tier"},
        {"stage": "Target-supported high v3", "count": m["target"]["high_target_supported_rows"], "note": "receiver target proxy support"},
        {"stage": "Top40 evidence review", "count": m["public"]["top_n_checked"], "note": "manual public evidence review"},
        {"stage": "Same-cancer axes", "count": m["case"]["same_cancer_rows"], "note": "safe case-study shortlist"},
    ]
    path = OUT_DIR / "global_overview_attrition.tsv"
    write_tsv(path, rows, ["stage", "count", "note"])
    tables["attrition"] = path

    rows = [
        {"method_or_union": "Marker baseline", "candidate_rows": m["pilot"]["candidates"]},
        {"method_or_union": "COMMOT", "candidate_rows": m["commot"]["candidate_rows"]},
        {"method_or_union": "LIANA", "candidate_rows": m["liana"]["candidate_rows"]},
        {"method_or_union": "Tri-method union", "candidate_rows": m["tri"]["union"]},
    ]
    path = OUT_DIR / "global_overview_method_candidate_counts.tsv"
    write_tsv(path, rows, ["method_or_union", "candidate_rows"])
    tables["method_counts"] = path

    rows = [
        {"support_count": "1 method", "candidate_keys": m["tri"]["union"] - m["tri"]["at_least_two"]},
        {"support_count": "2 methods", "candidate_keys": m["tri"]["at_least_two"] - m["tri"]["all_three"]},
        {"support_count": "3 methods", "candidate_keys": m["tri"]["all_three"]},
    ]
    path = OUT_DIR / "global_overview_method_support.tsv"
    write_tsv(path, rows, ["support_count", "candidate_keys"])
    tables["method_support"] = path

    rows = [
        {"null_model": "Spatial", "pass_p_le_0_05": m["spatial_null"]["p_le_0_05"], "tested": m["spatial_null"]["candidate_rows"]},
        {"null_model": "Label", "pass_p_le_0_05": m["label_null"]["p_le_0_05"], "tested": m["label_null"]["candidate_rows"]},
        {"null_model": "Fake L-R", "pass_p_le_0_05": m["fake_null"]["p_le_0_05"], "tested": m["fake_null"]["candidate_rows"]},
        {"null_model": "All gates", "pass_p_le_0_05": m["diag"]["null_all_pass_rows"], "tested": m["diag"]["rows"]},
    ]
    path = OUT_DIR / "global_overview_null_pass.tsv"
    write_tsv(path, rows, ["null_model", "pass_p_le_0_05", "tested"])
    tables["null_pass"] = path

    rows = [
        {"tier": "High", "count": m["score"]["tiers"]["high_pilot_v3"]},
        {"tier": "Medium", "count": m["score"]["tiers"]["medium_pilot_v3"]},
        {"tier": "Low", "count": m["score"]["tiers"]["low_pilot_v3"]},
    ]
    path = OUT_DIR / "global_overview_v3_tiers.tsv"
    write_tsv(path, rows, ["tier", "count"])
    tables["v3_tiers"] = path

    rows = [{"pathway": k, "high_candidates": v} for k, v in sorted(m["diag"]["top_high_pathways"].items(), key=lambda x: (-x[1], x[0]))]
    path = OUT_DIR / "global_overview_high_pathways.tsv"
    write_tsv(path, rows, ["pathway", "high_candidates"])
    tables["high_pathways"] = path

    rows = [
        {"group": "All scored", "supported": m["target"]["target_supported_rows"], "unsupported": m["target"]["joined_score_rows"] - m["target"]["target_supported_rows"]},
        {"group": "High v3", "supported": m["target"]["high_target_supported_rows"], "unsupported": m["gap"]["unsupported_high_v3_rows"]},
    ]
    path = OUT_DIR / "global_overview_target_support.tsv"
    write_tsv(path, rows, ["group", "supported", "unsupported"])
    tables["target_support"] = path

    gap_counts = m["gap"]["reason_counts"]
    rows = [
        {"gap_reason": "Low percentile + enrichment", "count": gap_counts.get("low_percentile_and_low_enrichment", 0)},
        {"gap_reason": "Low enrichment", "count": gap_counts.get("low_receiver_enrichment", 0)},
        {"gap_reason": "Low percentile", "count": gap_counts.get("low_receiver_percentile", 0)},
    ]
    path = OUT_DIR / "global_overview_target_gap_reasons.tsv"
    write_tsv(path, rows, ["gap_reason", "count"])
    tables["target_gaps"] = path

    rows = [{"category": k.replace("_", " "), "count": v} for k, v in m["public"]["manual_review_category_counts"].items()]
    path = OUT_DIR / "global_overview_public_evidence.tsv"
    write_tsv(path, rows, ["category", "count"])
    tables["public_evidence"] = path

    rows = []
    for name, value in [
        ("minimal_input_samples", len(m["qc_rows"])),
        ("marker_baseline_candidates", m["pilot"]["candidates"]),
        ("tri_method_union", m["tri"]["union"]),
        ("all_three_methods", m["tri"]["all_three"]),
        ("all_null_gates", m["diag"]["null_all_pass_rows"]),
        ("high_v3", m["score"]["tiers"]["high_pilot_v3"]),
        ("target_supported_high_v3", m["target"]["high_target_supported_rows"]),
        ("top40_public_evidence_review", m["public"]["top_n_checked"]),
        ("same_cancer_axes", m["case"]["same_cancer_rows"]),
        ("same_sample_validation_upgrades", m["audit"]["evidence_upgrade_counts"].get("exact_sample_or_original_paper", 0)),
    ]:
        rows.append({"metric": name, "value": value})
    path = OUT_DIR / "global_overview_key_metrics.tsv"
    write_tsv(path, rows, ["metric", "value"])
    tables["key_metrics"] = path
    return tables


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=9, fontweight="bold")


def style_axis(ax) -> None:
    ax.tick_params(axis="both", labelsize=6, length=2)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    ax.set_axisbelow(True)


def draw_attrition(ax, df: pd.DataFrame) -> None:
    ax.axis("off")
    ax.set_title("Formal scoring is separated from contextual prioritization", loc="left", fontsize=8.2, pad=2)
    xs = np.linspace(0.04, 0.96, len(df))
    colors = [PALETTE["gray"], PALETTE["light_blue"], PALETTE["blue"], PALETTE["green"], PALETTE["yellow"], PALETTE["orange"]]
    for i, row in df.iterrows():
        x = xs[i]
        w = 0.135
        h = 0.46
        box = FancyBboxPatch((x - w / 2, 0.31), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                             facecolor=colors[i], edgecolor="white", linewidth=1.0, transform=ax.transAxes)
        ax.add_patch(box)
        ax.text(x, 0.61, f"{int(row['count']):,}", transform=ax.transAxes, ha="center", va="center", fontsize=11, fontweight="bold", color="white" if i in [2, 3] else "#222222")
        label = row["stage"].replace("Target-supported high v3", "Target-supported\nhigh v3").replace("Top40 evidence review", "Top40 evidence\nreview").replace("Marker baseline", "Marker\nbaseline").replace("Same-cancer axes", "Same-cancer\naxes")
        ax.text(x, 0.43, label, transform=ax.transAxes, ha="center", va="center", fontsize=6.2, color="white" if i in [2, 3] else "#222222")
        if i < len(df) - 1:
            ax.annotate("", xy=(xs[i + 1] - w / 2 - 0.012, 0.54), xytext=(x + w / 2 + 0.012, 0.54),
                        xycoords=ax.transAxes, textcoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color="#666666", linewidth=0.8))
    divider_x = (xs[2] + xs[3]) / 2
    ax.plot([divider_x, divider_x], [0.24, 0.84], transform=ax.transAxes, color="#777777", linewidth=0.7, linestyle="--")
    ax.text(0.04, 0.84, "Formal score and gates", transform=ax.transAxes, ha="left", va="bottom", fontsize=6.0, fontweight="bold", color="#444444")
    ax.text(divider_x + 0.02, 0.84, "Contextual review; tier unchanged", transform=ax.transAxes, ha="left", va="bottom", fontsize=6.0, fontweight="bold", color="#444444")
    ax.text(0.04, 0.13, "Counts are computational prioritization outputs, not experimental validation.", transform=ax.transAxes, ha="left", va="center", fontsize=6.3, color="#555555")


def draw_method_support(ax, method_df: pd.DataFrame, support_df: pd.DataFrame) -> None:
    ax.set_title("Method support is informative but incomplete", loc="left", fontsize=7.6)
    y = np.arange(len(support_df))
    bars = ax.barh(y, support_df["candidate_keys"], color=[PALETTE["gray"], PALETTE["light_blue"], PALETTE["blue"]], height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(support_df["support_count"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Candidate keys", fontsize=6)
    style_axis(ax)
    for bar, val in zip(bars, support_df["candidate_keys"]):
        ax.text(val + 70, bar.get_y() + bar.get_height() / 2, f"{int(val):,}", va="center", fontsize=6)
    ax.text(0.98, 0.94, "Union: 7,346; all three: 2,924", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color="#555555")


def draw_null_pass(ax, df: pd.DataFrame) -> None:
    ax.set_title("Null screens use explicitly distinct criteria", loc="left", fontsize=7.6)
    colors = [PALETTE["light_blue"], PALETTE["purple"], PALETTE["green"], PALETTE["blue"]]
    y = np.arange(len(df))
    bars = ax.barh(y, df["pass_p_le_0_05"], color=colors, height=0.58)
    ax.set_yticks(y)
    labels = ["All gates\n(all p <= 0.10)" if value == "All gates" else f"{value}\n(p <= 0.05)" for value in df["null_model"]]
    ax.set_yticklabels(labels, fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("Candidates passing the stated criterion", fontsize=6)
    style_axis(ax)
    for bar, val in zip(bars, df["pass_p_le_0_05"]):
        ax.text(val + 35, bar.get_y() + bar.get_height() / 2, f"{int(val):,}", va="center", fontsize=6)


def draw_tiers(ax, df: pd.DataFrame) -> None:
    ax.set_title("v3 score assigns operational pilot tiers", loc="left", fontsize=7.6)
    order = ["High", "Medium", "Low"]
    colors = [PALETTE["blue"], PALETTE["light_blue"], PALETTE["gray"]]
    x0 = 0
    total = float(df["count"].sum())
    for label, color in zip(order, colors):
        val = int(df.loc[df["tier"] == label, "count"].iloc[0])
        ax.barh([0], [val], left=x0, color=color, height=0.42)
        if val < 420:
            ax.text(x0 + val + 80, 0.18, f"{label}: {val:,}", ha="left", va="bottom", fontsize=5.8, color="#222222")
        else:
            ax.text(x0 + val / 2, 0, f"{label}\n{val:,}", ha="center", va="center", fontsize=6, color="#222222")
        x0 += val
    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel("Scored candidates", fontsize=6)
    style_axis(ax)


def draw_target(ax, support_df: pd.DataFrame, gap_df: pd.DataFrame) -> None:
    ax.set_title("Receiver target proxy contextualizes high-v3 candidates", loc="left", fontsize=7.6)
    groups = support_df["group"].tolist()
    y = np.arange(len(groups))
    ax.barh(y, support_df["supported"], color=PALETTE["green"], height=0.42, label="Supported")
    ax.barh(y, support_df["unsupported"], left=support_df["supported"], color=PALETTE["gray"], height=0.42, label="Not supported")
    ax.set_yticks(y)
    ax.set_yticklabels(groups, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Candidates", fontsize=6)
    style_axis(ax)
    for yy, supp, unsupp in zip(y, support_df["supported"], support_df["unsupported"]):
        if supp < 350:
            ax.text(supp + 70, yy - 0.12, f"{int(supp):,}", ha="left", va="center", fontsize=5.8, color=PALETTE["green"])
        else:
            ax.text(supp / 2, yy, f"{int(supp):,}", ha="center", va="center", fontsize=6, color="white")
        if unsupp < 350:
            ax.text(supp + unsupp + 70, yy + 0.12, f"{int(unsupp):,}", ha="left", va="center", fontsize=5.8, color="#555555")
        else:
            ax.text(supp + unsupp / 2, yy, f"{int(unsupp):,}", ha="center", va="center", fontsize=6, color="#222222")
    ax.legend(loc="lower right", fontsize=5.8, handlelength=0.8)


def draw_pathways(ax, df: pd.DataFrame) -> None:
    ax.set_title("High-priority v3 candidates span canonical TME pathways", loc="left", fontsize=7.6)
    plot_df = df.sort_values("high_candidates", ascending=True)
    ax.barh(plot_df["pathway"], plot_df["high_candidates"], color=PALETTE["blue"], height=0.55)
    ax.set_xlabel("High-priority candidates", fontsize=6)
    style_axis(ax)
    ax.tick_params(axis="y", labelsize=5.8)
    for ytick, val in enumerate(plot_df["high_candidates"]):
        ax.text(val + 1, ytick, str(int(val)), va="center", fontsize=5.8)


def draw_public_evidence(ax, df: pd.DataFrame) -> None:
    ax.set_title("Public evidence supports contextual case selection", loc="left", fontsize=7.6)
    order = ["same cancer axis", "other cancer axis", "broad pathway"]
    colors = [PALETTE["green"], PALETTE["yellow"], PALETTE["gray"]]
    vals = [int(df.loc[df["category"] == item, "count"].iloc[0]) for item in order]
    y = np.arange(len(order))
    bars = ax.barh(y, vals, color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(["Same-cancer axis", "Other-cancer axis", "Broad pathway"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Top40 reviewed candidates", fontsize=6)
    style_axis(ax)
    for bar, val in zip(bars, vals):
        ax.text(val + 0.5, bar.get_y() + bar.get_height() / 2, str(val), va="center", fontsize=6)
    ax.text(0.02, 0.02, "0/7 upgraded to same-sample validation", transform=ax.transAxes, fontsize=6, color="#555555")


def render(tables: dict[str, Path], m: dict) -> dict:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    attrition = pd.read_csv(tables["attrition"], sep="\t")
    method_counts = pd.read_csv(tables["method_counts"], sep="\t")
    method_support = pd.read_csv(tables["method_support"], sep="\t")
    null_pass = pd.read_csv(tables["null_pass"], sep="\t")
    tiers = pd.read_csv(tables["v3_tiers"], sep="\t")
    target_support = pd.read_csv(tables["target_support"], sep="\t")
    target_gaps = pd.read_csv(tables["target_gaps"], sep="\t")
    high_pathways = pd.read_csv(tables["high_pathways"], sep="\t")
    public_evidence = pd.read_csv(tables["public_evidence"], sep="\t")

    fig = plt.figure(figsize=(9.4, 7.3), constrained_layout=True)
    axes = fig.subplot_mosaic(
        [["A", "A", "A"], ["B", "C", "D"], ["E", "F", "G"]],
        gridspec_kw={"height_ratios": [1.0, 1.18, 1.38]},
    )
    draw_attrition(axes["A"], attrition)
    draw_method_support(axes["B"], method_counts, method_support)
    draw_null_pass(axes["C"], null_pass)
    draw_tiers(axes["D"], tiers)
    draw_target(axes["E"], target_support, target_gaps)
    draw_pathways(axes["F"], high_pathways)
    draw_public_evidence(axes["G"], public_evidence)
    for label in "ABCDEFG":
        add_panel_label(axes[label], label.lower())

    fig.suptitle("SpatialLR-Trust global Results overview", x=0.01, y=1.02, ha="left", fontsize=10, fontweight="bold")
    base = FIG_DIR / BASE_NAME
    outputs = {
        "svg": str(base.with_suffix(".svg").relative_to(PROJECT)),
        "pdf": str(base.with_suffix(".pdf").relative_to(PROJECT)),
        "tiff": str(base.with_suffix(".tiff").relative_to(PROJECT)),
        "png": str(base.with_suffix(".png").relative_to(PROJECT)),
    }
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "run_id": RUN_ID,
        "figure": "SpatialLR-Trust global Results overview",
        "core_conclusion": "SpatialLR-Trust applies formal operational scoring before separate contextual target and public-evidence review.",
        "archetype": "quantitative grid",
        "backend": "python/matplotlib",
        "input_scope": "Derived TSV and JSON outputs only; no raw expression matrix was read.",
        "source_tables": {name: str(path.relative_to(PROJECT)) for name, path in tables.items()},
        "outputs": outputs,
        "key_counts": {
            "minimal_input_samples": len(m["qc_rows"]),
            "marker_baseline_candidates": m["pilot"]["candidates"],
            "tri_method_union": m["tri"]["union"],
            "all_three_methods": m["tri"]["all_three"],
            "all_null_gates": m["diag"]["null_all_pass_rows"],
            "high_v3": m["score"]["tiers"]["high_pilot_v3"],
            "target_supported_high_v3": m["target"]["high_target_supported_rows"],
            "top40_public_evidence_review": m["public"]["top_n_checked"],
            "same_cancer_axes": m["case"]["same_cancer_rows"],
            "same_sample_validation_upgrades": m["audit"]["evidence_upgrade_counts"].get("exact_sample_or_original_paper", 0),
        },
        "interpretation_boundary": "Computational benchmark figure; evidence supports prioritization and filtering, not experimental validation.",
    }
    manifest_path = FIG_DIR / f"{BASE_NAME}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with (RUN_DIR / "outputs_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["artifact", "path"])
        for name, path in outputs.items():
            writer.writerow([name, path])
        writer.writerow(["figure_manifest", str(manifest_path.relative_to(PROJECT))])
        for name, path in tables.items():
            writer.writerow(["source_" + name, str(path.relative_to(PROJECT))])
    return manifest


def main() -> None:
    m = load_metrics()
    tables = build_source_tables(m)
    manifest = render(tables, m)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
