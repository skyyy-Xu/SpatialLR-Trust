#!/usr/bin/env python3
"""Render a plot-ready MIF-CD74 case-study figure from derived TSV inputs."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "20260709_2305_mif-cd74-case-figure")
CASE_DIR = PROJECT / "results/task_e_case_studies/GSE292299_GSM8855712_MIF_CD74"
FIG_DIR = PROJECT / "figures/task_e_case_studies/GSE292299_GSM8855712_MIF_CD74"
RUN_DIR = PROJECT / "runs" / RUN_ID
TF48_DIR = PROJECT / "results/task_f_reproducibility/tf48_stage4_revision_analyses" / (
    "20260810_1631_tf48-stage4-revision-analyses"
)
BASE_NAME = "gsm8855712_mif_cd74_case_study_panel"

PANEL_LETTERS = ["a", "b", "c", "d", "e", "f"]
COMPARTMENT_COLORS = {
    "tumor_like": "#4C78A8",
    "immune": "#59A14F",
    "stroma_fibroblast": "#E15759",
    "endothelial": "#76B7B2",
    "ambiguous_or_low_signal": "#B9B9B9",
}
COMPARTMENT_LABELS = {
    "tumor_like": "Tumor-like",
    "immune": "Immune",
    "stroma_fibroblast": "Stroma/fibroblast",
    "endothelial": "Endothelial",
    "ambiguous_or_low_signal": "Ambiguous/low",
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


def read_evidence(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return dict(csv.reader(handle, delimiter="\t"))


def prep_spots(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["x"] = pd.to_numeric(out["pxl_col_in_fullres"])
    out["y"] = -pd.to_numeric(out["pxl_row_in_fullres"])
    for col in ["MIF", "CD74", "mif_target_mean", "mif_cd74_min_expr"]:
        out[col] = pd.to_numeric(out[col])
        out[f"log1p_{col}"] = np.log1p(out[col])
    out["mif_cd74_both_positive"] = pd.to_numeric(out["mif_cd74_both_positive"]).astype(int)
    return out


def add_panel_label(ax, letter: str) -> None:
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, ha="right", va="bottom", fontweight="bold", fontsize=9)


def style_spatial_axis(ax, title: str) -> None:
    ax.set_title(title, loc="left", pad=2, fontsize=8)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def scatter_continuous(fig, ax, df, value_col: str, title: str, cmap: str, cbar_label: str) -> None:
    order = np.argsort(df[value_col].to_numpy())
    vals = df[value_col].to_numpy()[order]
    vmax = np.percentile(vals, 99.5) if len(vals) else 1.0
    sc = ax.scatter(
        df["x"].to_numpy()[order],
        df["y"].to_numpy()[order],
        c=vals,
        s=1.8,
        cmap=cmap,
        vmin=0,
        vmax=max(vmax, 1e-9),
        linewidths=0,
        rasterized=True,
    )
    style_spatial_axis(ax, title)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.01)
    cbar.set_label(cbar_label, fontsize=6)
    cbar.ax.tick_params(labelsize=6, length=2)


def draw_compartment(ax, df, comp_summary) -> None:
    for comp, color in COMPARTMENT_COLORS.items():
        sub = df[df["spot_compartment"] == comp]
        if sub.empty:
            continue
        ax.scatter(sub["x"], sub["y"], s=1.8, c=color, linewidths=0, alpha=0.95, rasterized=True)
    style_spatial_axis(ax, "Coarse compartment map")


def draw_copositive(ax, df) -> None:
    ax.scatter(df["x"], df["y"], s=1.4, c="#D4D4D4", linewidths=0, alpha=0.65, rasterized=True)
    hit = df[df["mif_cd74_both_positive"] == 1]
    ax.scatter(hit["x"], hit["y"], s=2.0, c="#C44E52", linewidths=0, alpha=0.9, rasterized=True)
    style_spatial_axis(ax, "MIF/CD74 co-positive spots")
    ax.text(0.02, 0.04, f"co-positive spots: {len(hit):,}", transform=ax.transAxes, ha="left", va="bottom", fontsize=6,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5))


def draw_summary(ax, comp_summary, evidence, target_summary) -> None:
    ax.axis("off")
    inner = ax.inset_axes([0.02, 0.52, 0.96, 0.43])
    plot_summary = comp_summary.copy()
    plot_summary["mean_mif_target"] = pd.to_numeric(plot_summary["mean_mif_target"])
    colors = [COMPARTMENT_COLORS[c] for c in plot_summary["spot_compartment"]]
    y = np.arange(len(plot_summary))
    inner.barh(y, plot_summary["mean_mif_target"], color=colors, height=0.62)
    inner.set_yticks(y)
    inner.set_yticklabels([COMPARTMENT_LABELS[c] for c in plot_summary["spot_compartment"]], fontsize=6)
    inner.invert_yaxis()
    inner.set_title("Original MIF-panel mean (CD74-inclusive)", loc="left", fontsize=6.3, pad=1)
    inner.tick_params(axis="x", labelsize=6, length=2)
    inner.tick_params(axis="y", length=0)
    inner.spines["left"].set_visible(False)
    for yy, val in zip(y, plot_summary["mean_mif_target"]):
        inner.text(val + 0.03, yy, f"{val:.2f}", va="center", fontsize=6)
    inner.set_xlim(0, max(plot_summary["mean_mif_target"].max() * 1.24, 0.1))

    text_lines = [
        "Cautionary audit (not validation)",
        "GSM8855712_NSCLC_P10",
        "tumor-like -> tumor-like MIF-CD74",
        "",
        f"SpatialLR-Trust v3 score: {float(evidence['v3_score']):.3f}",
        f"MIF/CD74 both-positive tumor-like spots: {float(evidence['tumor_like_fraction_mif_cd74_both_positive']) * 100:.1f}%",
        "Original target panel includes CD74;",
        "the original enrichment is not independent.",
        f"MIF support after exclusion: {target_summary['mif_original_supported_rows']} -> {target_summary['mif_excluded_supported_rows']} / {target_summary['mif_rows']} rows",
        "",
        "Dominant spot compartment; directionality and",
        "same-sample experimental validation remain unresolved.",
    ]
    ax.text(0.02, 0.37, "\n".join(text_lines), transform=ax.transAxes, ha="left", va="top", fontsize=6.5, linespacing=1.22)


def render() -> dict[str, str]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    spots = prep_spots(pd.read_csv(CASE_DIR / "gsm8855712_mif_cd74_spot_plot_data.tsv", sep="\t"))
    comp_summary = pd.read_csv(CASE_DIR / "gsm8855712_mif_cd74_compartment_summary.tsv", sep="\t")
    evidence = read_evidence(CASE_DIR / "gsm8855712_mif_cd74_evidence_summary.tsv")
    target_summary = read_evidence(TF48_DIR / "target_candidate_excluded_summary.tsv")

    fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=True)
    mosaic = [["a", "b", "e"], ["c", "d", "f"]]
    axes = fig.subplot_mosaic(mosaic, width_ratios=[1, 1, 1.15], height_ratios=[1, 1])

    draw_compartment(axes["a"], spots, comp_summary)
    scatter_continuous(fig, axes["b"], spots, "log1p_MIF", "MIF expression", "magma", "log1p counts")
    scatter_continuous(fig, axes["c"], spots, "log1p_CD74", "CD74 expression", "viridis", "log1p counts")
    scatter_continuous(fig, axes["d"], spots, "log1p_mif_target_mean", "Original MIF-panel mean", "cividis", "log1p mean")
    draw_copositive(axes["e"], spots)
    draw_summary(axes["f"], comp_summary, evidence, target_summary)

    for letter, key in zip(PANEL_LETTERS, ["a", "b", "c", "d", "e", "f"]):
        add_panel_label(axes[key], letter)

    fig.suptitle("SpatialLR-Trust cautionary audit: MIF-CD74 in coarse NSCLC spot data", x=0.01, y=1.02,
                 ha="left", fontsize=9, fontweight="bold")

    base = FIG_DIR / BASE_NAME
    outputs = {
        "svg": str(base.with_suffix(".svg").relative_to(PROJECT)),
        "pdf": str(base.with_suffix(".pdf").relative_to(PROJECT)),
        "tiff": str(base.with_suffix(".tiff").relative_to(PROJECT)),
        "png": str(base.with_suffix(".png").relative_to(PROJECT)),
    }
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    svg_path = base.with_suffix(".svg")
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "run_id": RUN_ID,
        "figure": "MIF-CD74 cautionary audit panel",
        "core_conclusion": "A high-scoring MIF-CD74 row remains non-resolving in a dominant spot compartment, and its original CD74-inclusive target proxy is non-independent.",
        "archetype": "image plate + quant",
        "backend": "python/matplotlib",
        "input_tables": [
            str((CASE_DIR / "gsm8855712_mif_cd74_spot_plot_data.tsv").relative_to(PROJECT)),
            str((CASE_DIR / "gsm8855712_mif_cd74_compartment_summary.tsv").relative_to(PROJECT)),
            str((CASE_DIR / "gsm8855712_mif_cd74_evidence_summary.tsv").relative_to(PROJECT)),
        ],
        "outputs": outputs,
        "n_spots": int(len(spots)),
        "tumor_like_spots": int(comp_summary.loc[comp_summary["spot_compartment"] == "tumor_like", "n_spots"].iloc[0]),
        "original_receiver_target_enrichment_non_independent": evidence["receiver_target_enrichment_from_task_e"],
        "tumor_like_mif_cd74_both_positive_fraction": evidence["tumor_like_fraction_mif_cd74_both_positive"],
        "mif_rows": int(target_summary["mif_rows"]),
        "mif_original_supported_rows": int(target_summary["mif_original_supported_rows"]),
        "mif_excluded_supported_rows": int(target_summary["mif_excluded_supported_rows"]),
        "interpretation_boundary": "Cautionary dominant-compartment audit; original target proxy is non-independent; no directionality or same-sample experimental validation.",
    }
    (FIG_DIR / f"{BASE_NAME}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    with (RUN_DIR / "outputs_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["artifact", "path"])
        for key, path in outputs.items():
            writer.writerow([key, path])
        writer.writerow(["manifest", str((FIG_DIR / f"{BASE_NAME}_manifest.json").relative_to(PROJECT))])
    return manifest


def main() -> None:
    print(json.dumps(render(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
