#!/usr/bin/env python3
"""Plot the manuscript-facing F5 full benchmark evidence summary."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


METHOD_ORDER = ["commot", "liana", "cellchat", "spatiallr_trust"]
METHOD_LABELS = {
    "commot": "COMMOT",
    "liana": "LIANA",
    "cellchat": "CellChat",
    "spatiallr_trust": "SpatialLR-Trust",
}
COLORS = {
    "commot": "#3B6FB6",
    "liana": "#8064A2",
    "cellchat": "#D08A32",
    "spatiallr_trust": "#237A70",
}


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def annotate_heatmap(ax: plt.Axes, values: np.ndarray) -> None:
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            color = "white" if value >= 0.62 else "#252525"
            ax.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=6,
                color=color,
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    source = Path(args.score_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    bootstrap = pd.read_csv(source / "section_bootstrap_metrics.tsv", sep="\t")
    hard = pd.read_csv(
        source / "hard_negative_false_positive_rates.tsv", sep="\t"
    )
    monotonic = pd.read_csv(
        source / "perturbation_monotonicity_summary.tsv", sep="\t"
    )
    stability = pd.read_csv(
        source / "technical_replicate_rank_stability_summary.tsv", sep="\t"
    )

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(7.2, 6.7), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        2,
        left=0.095,
        right=0.985,
        bottom=0.09,
        top=0.965,
        width_ratios=[1.12, 1.0],
        height_ratios=[0.90, 1.15],
        wspace=0.34,
        hspace=0.42,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    lower_right = grid[1, 1].subgridspec(2, 1, hspace=0.82)
    ax_c = fig.add_subplot(lower_right[0, 0])
    ax_d = fig.add_subplot(lower_right[1, 0])

    primary = bootstrap[
        bootstrap["metric_universe"].eq("common_three_method_primary")
    ].copy()
    offsets = {"auroc": 0.13, "average_precision": -0.13}
    markers = {"auroc": "o", "average_precision": "s"}
    for y_index, scorer in enumerate(METHOD_ORDER):
        for metric_name in ["auroc", "average_precision"]:
            row = primary[
                primary["scorer"].eq(scorer)
                & primary["metric"].eq(metric_name)
            ].iloc[0]
            y = len(METHOD_ORDER) - 1 - y_index + offsets[metric_name]
            ax_a.errorbar(
                row["point_estimate"],
                y,
                xerr=np.array(
                    [
                        [row["point_estimate"] - row["ci_lower_0_025"]],
                        [row["ci_upper_0_975"] - row["point_estimate"]],
                    ]
                ),
                fmt=markers[metric_name],
                color=COLORS[scorer],
                markeredgecolor="white",
                markeredgewidth=0.4,
                markersize=5.2,
                elinewidth=1.1,
                capsize=2,
                zorder=3,
            )
    ax_a.axvline(0.5, color="#A0A0A0", linestyle="--", linewidth=0.8)
    ax_a.set_xlim(0.44, 0.90)
    ax_a.set_ylim(-0.55, 3.55)
    ax_a.set_yticks(range(4))
    ax_a.set_yticklabels(
        [METHOD_LABELS[item] for item in reversed(METHOD_ORDER)]
    )
    ax_a.set_xlabel(
        "Discrimination on the common three-method universe "
        "(section-cluster 95% bootstrap interval)"
    )
    ax_a.set_title(
        "Frozen score shows a modest, uncertain gain on the common universe",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=7,
    )
    ax_a.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    ax_a.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="#555555",
                linestyle="none",
                label="AUROC",
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                color="#555555",
                linestyle="none",
                label="Average precision",
            ),
        ],
        loc="lower right",
        ncol=2,
        fontsize=6.5,
    )
    panel_label(ax_a, "a")

    selection_order = [
        "commot_detected",
        "liana_detected",
        "cellchat_detected",
        "trust_nonlow",
        "trust_high",
    ]
    selection_labels = [
        "COMMOT detected",
        "LIANA detected",
        "CellChat detected",
        "Trust non-low",
        "Trust high",
    ]
    hard_order = [
        "hard_negative_ligand_only",
        "hard_negative_receptor_only",
        "hard_negative_reporter_absent",
        "hard_negative_spatial_separated",
    ]
    hard_labels = [
        "Ligand\nonly",
        "Receptor\nonly",
        "Reporter\nabsent",
        "Spatially\nseparated",
    ]
    hard_matrix = (
        hard.pivot(
            index="selection",
            columns="condition_id",
            values="false_positive_rate",
        )
        .loc[selection_order, hard_order]
        .to_numpy()
    )
    image_b = ax_b.imshow(
        hard_matrix, vmin=0, vmax=1, cmap="YlOrBr", aspect="auto"
    )
    annotate_heatmap(ax_b, hard_matrix)
    ax_b.set_xticks(range(4), hard_labels)
    ax_b.set_yticks(range(5), selection_labels)
    ax_b.set_title(
        "Specific hard negatives remain unresolved",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=7,
    )
    ax_b.set_xlabel("False-positive rate (n = 60 rows per family)")
    ax_b.tick_params(length=0)
    colorbar_b = fig.colorbar(image_b, ax=ax_b, fraction=0.035, pad=0.025)
    colorbar_b.set_label("FPR", fontsize=6)
    colorbar_b.ax.tick_params(labelsize=6)
    panel_label(ax_b, "b")

    family_order = [
        "coordinate_corruption",
        "label_corruption",
        "expression_corruption",
        "spike_strength",
    ]
    mono_matrix = (
        monotonic.pivot(
            index="scorer",
            columns="perturbation_family",
            values="fully_monotonic_fraction",
        )
        .loc[METHOD_ORDER, family_order]
        .to_numpy()
    )
    image_c = ax_c.imshow(
        mono_matrix, vmin=0, vmax=1, cmap="Blues", aspect="auto"
    )
    annotate_heatmap(ax_c, mono_matrix)
    ax_c.set_xticks(
        range(4), ["Coord.", "Label", "Expr.", "Spike"], rotation=0
    )
    ax_c.set_yticks(
        range(4), [METHOD_LABELS[item] for item in METHOD_ORDER]
    )
    ax_c.set_title(
        "Full-series monotonic response",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=6,
    )
    ax_c.tick_params(length=0)
    colorbar_c = fig.colorbar(image_c, ax=ax_c, fraction=0.045, pad=0.025)
    colorbar_c.set_label("Fraction", fontsize=6)
    colorbar_c.ax.tick_params(labelsize=6)
    panel_label(ax_c, "c")

    x_positions = np.arange(len(METHOD_ORDER))
    sample_ids = sorted(stability["sample_id"].unique())
    sample_offsets = np.linspace(-0.14, 0.14, len(sample_ids))
    for sample_offset, sample_id in zip(sample_offsets, sample_ids):
        section = stability[stability["sample_id"].eq(sample_id)].set_index(
            "scorer"
        )
        values = [section.loc[item, "mean_spearman"] for item in METHOD_ORDER]
        ax_d.scatter(
            x_positions + sample_offset,
            values,
            s=20,
            facecolors=[COLORS[item] for item in METHOD_ORDER],
            edgecolors="white",
            linewidths=0.4,
            alpha=0.9,
            zorder=3,
        )
    means = (
        stability.groupby("scorer")["mean_spearman"].mean().reindex(METHOD_ORDER)
    )
    ax_d.plot(
        x_positions,
        means,
        marker="_",
        markersize=13,
        linewidth=0,
        color="#202020",
        markeredgewidth=1.4,
        zorder=4,
    )
    ax_d.set_xticks(
        x_positions,
        ["COMMOT", "LIANA", "CellChat", "Trust"],
        rotation=20,
        ha="right",
    )
    ax_d.set_ylim(0.4, 1.02)
    ax_d.set_ylabel("Mean pairwise Spearman")
    ax_d.set_title(
        "Trust score is less stable across technical replicates",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=6,
    )
    ax_d.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    panel_label(ax_d, "d")

    base = output / "figure_f5_full_benchmark"
    svg_path = base.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(base.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    summary = {
        "status": "PASS_F5_FULL_BENCHMARK_FIGURE_BUILT",
        "run_id": args.run_id,
        "backend": "Python matplotlib",
        "archetype": "quantitative grid with hero comparison panel",
        "core_conclusion": (
            "The frozen score has a modest common-universe discrimination "
            "gain but retains hard-negative and stability limitations."
        ),
        "panels": {
            "a": "Common-universe AUROC and average precision with section bootstrap",
            "b": "Hard-negative family false-positive rates",
            "c": "Perturbation-series monotonicity",
            "d": "Technical-replicate rank stability",
        },
        "exports": ["svg", "pdf", "tiff", "png"],
    }
    (output / "figure_build_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
