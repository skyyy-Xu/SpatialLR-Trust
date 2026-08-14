#!/usr/bin/env python3
"""Generate one frozen 19-condition F5 section/replicate batch."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from generate_task_f_f5_semisim_smoke import (
    CELLS_PER_REGION,
    add_signal,
    median_cross_distance,
    mix_coordinates,
    mix_labels,
    nearest_region,
    rng_for,
    separated_region,
    write_condition,
)


EXPECTED_CONDITIONS = {
    "unmodified_reference",
    "spike_low",
    "spike_medium",
    "spike_high",
    "coordinate_mix_10",
    "coordinate_mix_25",
    "coordinate_mix_50",
    "coordinate_mix_100",
    "label_mix_10",
    "label_mix_25",
    "label_mix_50",
    "expression_mix_10",
    "expression_mix_25",
    "expression_mix_50",
    "hard_negative_spatial_separated",
    "hard_negative_ligand_only",
    "hard_negative_receptor_only",
    "hard_negative_reporter_absent",
    "fixed_fake_lr",
}


def mix_expression_with_audit(
    matrix: np.ndarray, fraction: float, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Permute selected cell profiles and return rows changed by that step."""
    if fraction <= 0:
        return matrix, np.zeros(len(matrix), dtype=bool)
    count = int(round(len(matrix) * fraction))
    selected = np.sort(rng.choice(len(matrix), size=count, replace=False))
    mixed = matrix.copy()
    mixed[selected] = matrix[rng.permutation(selected)]
    changed = np.any(mixed != matrix, axis=1)
    return mixed, changed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--pair-contract", required=True)
    parser.add_argument("--full-manifest", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--generated-root", required=True)
    parser.add_argument("--generated-path-prefix", required=True)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    generated_root = Path(args.generated_root).resolve()
    result_dir = Path(args.result_dir).resolve()
    generated_root.mkdir(parents=True, exist_ok=False)
    result_dir.mkdir(parents=True, exist_ok=True)
    data = ad.read_h5ad(input_dir / "shared_input.h5ad")
    pairs = pd.read_csv(args.pair_contract, sep="\t")
    full_manifest = pd.read_csv(args.full_manifest, sep="\t")
    manifest = full_manifest[
        full_manifest["sample_id"].eq(args.sample_id)
        & full_manifest["replicate"].eq(args.replicate)
    ].copy()
    if data.shape != (5250, 280) or len(pairs) != 4 or len(manifest) != 19:
        raise ValueError("F5 full-batch source cardinality mismatch")
    if set(manifest["condition_id"]) != EXPECTED_CONDITIONS:
        raise ValueError("F5 full-batch condition set mismatch")
    if set(data.obs["sample_id"].astype(str)) != {args.sample_id}:
        raise ValueError("F5 full-batch sample identity mismatch")
    if manifest["seed"].nunique() != 1:
        raise ValueError("F5 full-batch seed is not unique")
    data.uns["f5_semisim_sample_id"] = args.sample_id
    data.uns["f5_semisim_replicate"] = args.replicate

    base_matrix = data.X.toarray() if sparse.issparse(data.X) else np.asarray(data.X)
    base_matrix = base_matrix.astype(np.float64)
    base_labels = data.obs["harmonized_compartment"].astype(str).to_numpy()
    base_xy = np.asarray(data.obsm["spatial"], dtype=float)
    genes = data.var_names.astype(str).str.upper().tolist()
    gene_index = {gene: index for index, gene in enumerate(genes)}
    base_sd = base_matrix.std(axis=0, ddof=0)
    seed = int(manifest["seed"].iloc[0])

    adjacent: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    separated: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for row in pairs.itertuples(index=False):
        adjacent[row.axis_id] = nearest_region(
            base_xy,
            base_labels,
            row.sender_compartment,
            row.receiver_compartment,
            seed,
            row.axis_id,
        )
        separated[row.axis_id] = separated_region(
            base_xy,
            base_labels,
            row.sender_compartment,
            row.receiver_compartment,
        )

    truth_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    file_rows: list[dict[str, object]] = []
    for condition in manifest.sort_values("condition_id").itertuples(index=False):
        matrix = base_matrix.copy()
        labels = base_labels.copy()
        xy = base_xy.copy()
        condition_rng = rng_for(seed, condition.condition_id, "condition")
        total_added = 0.0
        for pair in pairs.itertuples(index=False):
            sender_cells, receiver_cells = adjacent[pair.axis_id]
            ligand = pair.ligand
            receptor = pair.receptor
            reporter = pair.synthetic_reporter
            mode = "none"
            if condition.condition_family in {
                "positive_spike",
                "coordinate_mix",
                "label_mix",
                "expression_mix",
            }:
                mode = "full"
            elif condition.condition_id == "hard_negative_spatial_separated":
                sender_cells, receiver_cells = separated[pair.axis_id]
                mode = "full"
            elif condition.condition_id == "hard_negative_ligand_only":
                mode = "ligand_only"
            elif condition.condition_id in {
                "hard_negative_receptor_only",
                "hard_negative_reporter_absent",
            }:
                mode = "full"
            elif condition.condition_id == "fixed_fake_lr":
                ligand = pair.fake_ligand
                receptor = pair.fake_receptor
                mode = "full"

            additions = {
                "ligand_addition": 0.0,
                "receptor_addition": 0.0,
                "reporter_addition": 0.0,
            }
            if mode != "none":
                additions = add_signal(
                    matrix,
                    base_sd,
                    gene_index,
                    ligand,
                    receptor,
                    reporter,
                    sender_cells,
                    receiver_cells,
                    float(condition.spike_strength_sd),
                    mode,
                )
                if condition.condition_id == "hard_negative_receptor_only":
                    matrix[sender_cells, gene_index[ligand]] -= additions[
                        "ligand_addition"
                    ]
                    additions["ligand_addition"] = 0.0
                if condition.condition_id == "hard_negative_reporter_absent":
                    matrix[receiver_cells, gene_index[reporter]] -= additions[
                        "reporter_addition"
                    ]
                    additions["reporter_addition"] = 0.0
                total_added += (
                    additions["ligand_addition"] * len(sender_cells)
                    + additions["receptor_addition"] * len(receiver_cells)
                    + additions["reporter_addition"] * len(receiver_cells)
                )
            truth_rows.append(
                {
                    "sample_id": args.sample_id,
                    "replicate": args.replicate,
                    "condition_id": condition.condition_id,
                    "condition_family": condition.condition_family,
                    "truth_class": condition.truth_class,
                    "axis_id": pair.axis_id,
                    "sender_compartment": pair.sender_compartment,
                    "receiver_compartment": pair.receiver_compartment,
                    "ligand": ligand,
                    "receptor": receptor,
                    "synthetic_reporter": reporter,
                    "resource_axis": int(condition.condition_id != "fixed_fake_lr"),
                    "sender_cells": len(sender_cells),
                    "receiver_cells": len(receiver_cells),
                    "median_cross_distance_before_corruption": median_cross_distance(
                        base_xy, sender_cells, receiver_cells
                    ),
                    **additions,
                }
            )

        xy = mix_coordinates(
            xy, float(condition.coordinate_mix_fraction), condition_rng
        )
        labels = mix_labels(
            labels, float(condition.label_mix_fraction), condition_rng
        )
        matrix, expression_changed = mix_expression_with_audit(
            matrix, float(condition.expression_mix_fraction), condition_rng
        )
        output = generated_root / condition.condition_id
        hashes = write_condition(
            output,
            data,
            matrix,
            labels,
            xy,
            condition.condition_id,
            input_dir,
        )
        matrix_changed = np.any(matrix != base_matrix, axis=1)
        coordinate_changed = np.any(xy != base_xy, axis=1)
        label_changed = labels != base_labels
        qc_rows.append(
            {
                "sample_id": args.sample_id,
                "replicate": args.replicate,
                "condition_id": condition.condition_id,
                "condition_family": condition.condition_family,
                "cells": len(matrix),
                "genes": matrix.shape[1],
                "matrix_changed_cells": int(matrix_changed.sum()),
                "matrix_changed_fraction": float(matrix_changed.mean()),
                "expression_changed_cells": int(expression_changed.sum()),
                "expression_changed_fraction": float(expression_changed.mean()),
                "coordinate_changed_cells": int(coordinate_changed.sum()),
                "coordinate_changed_fraction": float(coordinate_changed.mean()),
                "label_changed_cells": int(label_changed.sum()),
                "label_changed_fraction": float(label_changed.mean()),
                "minimum_expression": float(matrix.min()),
                "total_expression_delta": float(matrix.sum() - base_matrix.sum()),
                "planned_signal_addition": total_added,
            }
        )
        for role, sha256 in hashes.items():
            file_rows.append(
                {
                    "sample_id": args.sample_id,
                    "replicate": args.replicate,
                    "condition_id": condition.condition_id,
                    "role": role,
                    "path": (
                        f"{args.generated_path_prefix.rstrip('/')}/"
                        f"{condition.condition_id}/{role}"
                    ),
                    "bytes": (output / role).stat().st_size,
                    "sha256": sha256,
                }
            )

    truth = pd.DataFrame(truth_rows).sort_values(["condition_id", "axis_id"])
    qc = pd.DataFrame(qc_rows).sort_values("condition_id")
    files = pd.DataFrame(file_rows).sort_values(["condition_id", "role"])
    truth.to_csv(result_dir / "truth_contract.tsv", sep="\t", index=False)
    qc.to_csv(result_dir / "condition_qc.tsv", sep="\t", index=False)
    files.to_csv(result_dir / "generated_input_manifest.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS_F5_SEMISIM_FULL_BATCH_INPUTS_BUILT",
        "run_id": args.run_id,
        "sample_id": args.sample_id,
        "replicate": args.replicate,
        "conditions": len(manifest),
        "condition_axis_truth_rows": len(truth),
        "generated_files": len(files),
        "cells_per_condition": data.n_obs,
        "genes_per_condition": data.n_vars,
        "positive_axes": len(pairs),
        "cells_per_injected_region": CELLS_PER_REGION,
        "seed": seed,
        "score_formula_modified": False,
        "boundary": "Generated input batch only; no method or performance result.",
    }
    (result_dir / "build_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
