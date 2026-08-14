#!/usr/bin/env python3
"""Run a minimal COMMOT baseline on selected SpatialLR-Trust samples.

This is the first external-method Task C bridge. It converts the project minimal
Visium inputs into AnnData, runs COMMOT on the curated pilot LR panel, and emits
standardized compartment-to-compartment candidate tables.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import anndata as ad
import commot as ct
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "manual")
COMPARTMENT_ORDER = [
    "tumor_like",
    "immune",
    "stroma_fibroblast",
    "endothelial",
    "ambiguous_or_low_signal",
]
METHOD = "commot_pilot_v1"


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace", newline="")
    return path.open("rt", errors="replace", newline="")


def read_table(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_barcodes(path: Path) -> List[str]:
    with open_text(path) as handle:
        return [line.rstrip("\n").split("\t", 1)[0] for line in handle if line.strip()]


def unique_names(names: Iterable[str]) -> List[str]:
    seen: Dict[str, int] = defaultdict(int)
    out: List[str] = []
    for raw in names:
        name = str(raw).upper()
        seen[name] += 1
        out.append(name if seen[name] == 1 else f"{name}__DUP{seen[name]}")
    return out


def read_features(path: Path) -> Tuple[List[str], Dict[str, str]]:
    var_names: List[str] = []
    gene_id_by_symbol: Dict[str, str] = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[0] if fields else ""
            symbol = fields[1] if len(fields) > 1 and fields[1] else gene_id
            symbol_up = symbol.upper()
            var_names.append(symbol_up)
            gene_id_by_symbol.setdefault(symbol_up, gene_id)
    return unique_names(var_names), gene_id_by_symbol


def read_spatial(path: Path) -> pd.DataFrame:
    rows = []
    with open_text(path) as handle:
        for idx, line in enumerate(handle):
            fields = line.rstrip("\n").split(",")
            if not fields or fields[0] == "":
                continue
            if idx == 0 and "barcode" in fields:
                continue
            if len(fields) < 4:
                continue
            rows.append(
                {
                    "barcode": fields[0],
                    "in_tissue": int(float(fields[1])),
                    "array_row": float(fields[2]),
                    "array_col": float(fields[3]),
                }
            )
    return pd.DataFrame(rows).set_index("barcode")


def read_annotations(path: Path) -> pd.DataFrame:
    with open_text(path) as handle:
        df = pd.read_csv(handle, sep="\t")
    return df.set_index("barcode")


def read_lr_panel(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    out = pd.DataFrame(
        {
            0: df["ligand"].str.upper(),
            1: df["receptor"].str.upper(),
            2: df["pathway"],
        }
    )
    return out.drop_duplicates().reset_index(drop=True)


def choose_samples(rows: List[Dict[str, str]], sample_ids: List[str], max_samples: int) -> List[Dict[str, str]]:
    if sample_ids:
        wanted = set(sample_ids)
        selected = [row for row in rows if row["sample_id"] in wanted]
        missing = sorted(wanted - {row["sample_id"] for row in selected})
        if missing:
            raise SystemExit(f"Missing sample_id(s) in manifest: {', '.join(missing)}")
        return selected
    return rows[:max_samples]


def build_adata(row: Dict[str, str], annot_row: Dict[str, str], lr_df: pd.DataFrame, min_total_counts: int) -> ad.AnnData:
    matrix_path = PROJECT / row["expression_path"]
    barcodes_path = PROJECT / row["barcodes_path"]
    features_path = PROJECT / row["features_path"]
    spatial_path = PROJECT / row["spatial_positions_path"]
    annotation_path = PROJECT / annot_row["annotation_path"]

    barcodes = read_barcodes(barcodes_path)
    var_names, gene_id_by_symbol = read_features(features_path)
    X = mmread(str(matrix_path)).tocsr().transpose().tocsr()
    if X.shape[0] != len(barcodes):
        raise ValueError(f"Barcode count mismatch for {row['sample_id']}: matrix {X.shape[0]} vs barcodes {len(barcodes)}")
    if X.shape[1] != len(var_names):
        raise ValueError(f"Feature count mismatch for {row['sample_id']}: matrix {X.shape[1]} vs features {len(var_names)}")

    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    spatial = read_spatial(spatial_path)
    annot = read_annotations(annotation_path)
    obs = obs.join(spatial, how="left").join(annot, how="left")
    obs["spot_compartment"] = obs["spot_compartment"].fillna("unannotated")
    keep = (obs["in_tissue"].fillna(0).astype(int) == 1) & obs["spot_compartment"].isin(COMPARTMENT_ORDER)
    if min_total_counts > 0:
        totals = np.asarray(X.sum(axis=1)).ravel()
        keep = keep & (totals >= min_total_counts)
    keep_idx = np.where(keep.to_numpy())[0]
    if len(keep_idx) < 20:
        raise ValueError(f"Too few usable annotated spots for {row['sample_id']}: {len(keep_idx)}")

    X = X[keep_idx, :]
    obs = obs.iloc[keep_idx].copy()
    var = pd.DataFrame(index=pd.Index(var_names, name="gene_symbol"))
    var["gene_id"] = [gene_id_by_symbol.get(name.split("__DUP", 1)[0], "") for name in var_names]

    # COMMOT expects nonnegative expression and adata.obsm['spatial'].
    adata = ad.AnnData(X=X.astype(np.float64), obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = obs[["array_col", "array_row"]].to_numpy(dtype=float)
    adata.obs["spot_compartment"] = pd.Categorical(adata.obs["spot_compartment"].astype(str), categories=COMPARTMENT_ORDER)
    adata.uns["spatiallr_sample"] = {
        "dataset": row["dataset"],
        "sample_id": row["sample_id"],
        "cancer": row["cancer"],
        "platform": row["platform"],
        "genes_in_lr_panel": sorted(set(lr_df[0]).union(set(lr_df[1])).intersection(set(adata.var_names))),
    }
    return adata


def summarize_pair_matrix(adata: ad.AnnData, obsp_key: str) -> List[Tuple[str, str, float, int, int]]:
    mat = adata.obsp[obsp_key]
    if not sparse.issparse(mat):
        mat = sparse.csr_matrix(mat)
    labels = adata.obs["spot_compartment"].astype(str).to_numpy()
    out = []
    for sender in COMPARTMENT_ORDER:
        sender_idx = np.where(labels == sender)[0]
        if len(sender_idx) == 0:
            continue
        for receiver in COMPARTMENT_ORDER:
            receiver_idx = np.where(labels == receiver)[0]
            if len(receiver_idx) == 0:
                continue
            block = mat[sender_idx, :][:, receiver_idx]
            mean_score = float(block.mean()) if block.shape[0] and block.shape[1] else 0.0
            if mean_score > 0:
                out.append((sender, receiver, mean_score, int(len(sender_idx)), int(len(receiver_idx))))
    return out


def run_sample(row: Dict[str, str], annot_row: Dict[str, str], lr_df: pd.DataFrame, args: argparse.Namespace):
    adata = build_adata(row, annot_row, lr_df, args.min_total_counts)
    available = set(adata.var_names)
    sample_lr = lr_df[lr_df[0].isin(available) & lr_df[1].isin(available)].copy().reset_index(drop=True)
    if sample_lr.empty:
        raise ValueError(f"No pilot LR pairs available in expression matrix for {row['sample_id']}")

    ct.tl.spatial_communication(
        adata,
        database_name="pilotLR",
        df_ligrec=sample_lr,
        pathway_sum=True,
        heteromeric=False,
        dis_thr=args.distance_threshold,
        cot_nitermax=args.cot_nitermax,
    )

    candidate_rows = []
    for _, lr in sample_lr.iterrows():
        ligand, receptor, pathway = str(lr[0]), str(lr[1]), str(lr[2])
        obsp_key = f"commot-pilotLR-{ligand}-{receptor}"
        if obsp_key not in adata.obsp:
            continue
        for sender, receiver, score, n_sender, n_receiver in summarize_pair_matrix(adata, obsp_key):
            candidate_rows.append(
                {
                    "method": METHOD,
                    "dataset": row["dataset"],
                    "sample_id": row["sample_id"],
                    "cancer": row["cancer"],
                    "platform": row["platform"],
                    "sender_compartment": sender,
                    "receiver_compartment": receiver,
                    "ligand": ligand,
                    "receptor": receptor,
                    "pathway": pathway,
                    "commot_mean_score": score,
                    "sender_spots": n_sender,
                    "receiver_spots": n_receiver,
                    "distance_threshold_array_units": args.distance_threshold,
                }
            )

    positive = [r for r in candidate_rows if r["commot_mean_score"] > 0]
    score_values = np.array([r["commot_mean_score"] for r in positive], dtype=float) if positive else np.array([], dtype=float)
    sample_summary = {
        "method": METHOD,
        "dataset": row["dataset"],
        "sample_id": row["sample_id"],
        "cancer": row["cancer"],
        "platform": row["platform"],
        "spots_used": int(adata.n_obs),
        "genes_used": int(adata.n_vars),
        "lr_pairs_requested": int(lr_df.shape[0]),
        "lr_pairs_available": int(sample_lr.shape[0]),
        "candidate_rows": int(len(candidate_rows)),
        "positive_candidate_rows": int(len(positive)),
        "mean_positive_score": float(score_values.mean()) if score_values.size else 0.0,
        "max_positive_score": float(score_values.max()) if score_values.size else 0.0,
        "distance_threshold_array_units": args.distance_threshold,
        "cot_nitermax": args.cot_nitermax,
    }
    return candidate_rows, sample_summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run COMMOT pilot baseline for SpatialLR-Trust minimal inputs.")
    parser.add_argument("--max-samples", type=int, default=1)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--distance-threshold", type=float, default=2.5)
    parser.add_argument("--cot-nitermax", type=int, default=2000)
    parser.add_argument("--min-total-counts", type=int, default=1)
    parser.add_argument("--output-prefix", default="commot_pilot")
    args = parser.parse_args()

    manifest_rows = read_table(PROJECT / "docs/minimal_input_manifest.tsv")
    annot_rows = {row["sample_id"]: row for row in read_table(PROJECT / "docs/spot_annotation_manifest.tsv")}
    lr_df = read_lr_panel(PROJECT / "docs/pilot_lr_panel.tsv")
    selected = choose_samples(manifest_rows, args.sample_id, args.max_samples)

    out_dir = PROJECT / "results/task_c_commot_baseline"
    run_dir = PROJECT / "runs" / RUN_ID
    out_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    all_candidates = []
    summaries = []
    failures = []
    for row in selected:
        try:
            candidates, summary = run_sample(row, annot_rows[row["sample_id"]], lr_df, args)
            all_candidates.extend(candidates)
            summaries.append(summary)
        except Exception as exc:  # keep pilot failures inspectable across samples
            failures.append(
                {
                    "dataset": row.get("dataset", ""),
                    "sample_id": row.get("sample_id", ""),
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )

    prefix = args.output_prefix
    cand_path = out_dir / f"{prefix}_candidates.tsv"
    summ_path = out_dir / f"{prefix}_sample_summary.tsv"
    fail_path = out_dir / f"{prefix}_failures.tsv"
    json_path = out_dir / f"{prefix}_summary.json"

    pd.DataFrame(all_candidates).to_csv(cand_path, sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(summ_path, sep="\t", index=False)
    pd.DataFrame(failures).to_csv(fail_path, sep="\t", index=False)
    summary_json = {
        "run_id": RUN_ID,
        "method": METHOD,
        "samples_requested": len(selected),
        "samples_succeeded": len(summaries),
        "samples_failed": len(failures),
        "candidate_rows": len(all_candidates),
        "outputs": {
            "candidates": rel(cand_path),
            "sample_summary": rel(summ_path),
            "failures": rel(fail_path),
        },
        "parameters": {
            "distance_threshold_array_units": args.distance_threshold,
            "cot_nitermax": args.cot_nitermax,
            "max_samples": args.max_samples,
            "sample_id": args.sample_id,
            "output_prefix": args.output_prefix,
        },
        "note": "COMMOT pilot baseline using coarse computational spot compartments; not a final benchmark result.",
    }
    json_path.write_text(json.dumps(summary_json, indent=2, sort_keys=True) + "\n")

    manifest_path = run_dir / "outputs_manifest.tsv"
    with manifest_path.open("w") as handle:
        handle.write("path\ttype\tdescription\n")
        handle.write(f"{rel(cand_path)}\ttsv\tStandardized COMMOT pilot compartment LR candidates\n")
        handle.write(f"{rel(summ_path)}\ttsv\tSample-level COMMOT pilot summary\n")
        handle.write(f"{rel(fail_path)}\ttsv\tCOMMOT pilot failures, if any\n")
        handle.write(f"{rel(json_path)}\tjson\tCOMMOT pilot run summary\n")
        handle.write(f"{rel(manifest_path)}\ttsv\tRun output manifest\n")
    print(json.dumps(summary_json, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
