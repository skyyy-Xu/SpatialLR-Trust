#!/usr/bin/env python3
"""Run a LIANA CellPhoneDB-style pilot baseline on SpatialLR-Trust samples."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "manual")
METHOD = "liana_cellphonedb_pilot"
COMPARTMENT_ORDER = [
    "tumor_like",
    "immune",
    "stroma_fibroblast",
    "endothelial",
    "ambiguous_or_low_signal",
]


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
    seen: Dict[str, int] = {}
    out: List[str] = []
    for raw in names:
        name = str(raw).upper()
        seen[name] = seen.get(name, 0) + 1
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
    return pd.DataFrame({"ligand": df["ligand"].str.upper(), "receptor": df["receptor"].str.upper(), "pathway": df["pathway"]}).drop_duplicates()


def choose_samples(rows: List[Dict[str, str]], sample_ids: List[str], max_samples: int) -> List[Dict[str, str]]:
    if sample_ids:
        wanted = set(sample_ids)
        selected = [row for row in rows if row["sample_id"] in wanted]
        missing = sorted(wanted - {row["sample_id"] for row in selected})
        if missing:
            raise SystemExit(f"Missing sample_id(s) in manifest: {', '.join(missing)}")
        return selected
    return rows[:max_samples]


def build_adata(row: Dict[str, str], annot_row: Dict[str, str], lr_df: pd.DataFrame, min_total_counts: int, min_cells: int) -> ad.AnnData:
    matrix_path = PROJECT / row["expression_path"]
    barcodes_path = PROJECT / row["barcodes_path"]
    features_path = PROJECT / row["features_path"]
    spatial_path = PROJECT / row["spatial_positions_path"]
    annotation_path = PROJECT / annot_row["annotation_path"]

    barcodes = read_barcodes(barcodes_path)
    var_names, gene_id_by_symbol = read_features(features_path)
    X = mmread(str(matrix_path)).tocsr().transpose().tocsr()
    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    obs = obs.join(read_spatial(spatial_path), how="left").join(read_annotations(annotation_path), how="left")
    obs["spot_compartment"] = obs["spot_compartment"].fillna("unannotated")
    keep = (obs["in_tissue"].fillna(0).astype(int) == 1) & obs["spot_compartment"].isin(COMPARTMENT_ORDER)
    if min_total_counts > 0:
        keep = keep & (np.asarray(X.sum(axis=1)).ravel() >= min_total_counts)
    keep_idx = np.where(keep.to_numpy())[0]
    X = X[keep_idx, :]
    obs = obs.iloc[keep_idx].copy()

    counts = obs["spot_compartment"].value_counts()
    keep_compartments = set(counts[counts >= min_cells].index.astype(str))
    obs_keep = obs["spot_compartment"].isin(keep_compartments).to_numpy()
    if obs_keep.sum() < 20:
        raise ValueError(f"Too few spots after min_cells filter for {row['sample_id']}: {obs_keep.sum()}")
    X = X[obs_keep, :]
    obs = obs.iloc[np.where(obs_keep)[0]].copy()

    var = pd.DataFrame(index=pd.Index(var_names, name="gene_symbol"))
    var["gene_id"] = [gene_id_by_symbol.get(name.split("__DUP", 1)[0], "") for name in var_names]
    adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obs["spot_compartment"] = adata.obs["spot_compartment"].astype(str)
    adata.obsm["spatial"] = obs[["array_col", "array_row"]].to_numpy(dtype=float)

    available_lr_genes = set(lr_df["ligand"]).union(set(lr_df["receptor"])).intersection(set(adata.var_names))
    if not available_lr_genes:
        raise ValueError(f"No LR panel genes available for {row['sample_id']}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.uns["spatiallr_sample"] = {"dataset": row["dataset"], "sample_id": row["sample_id"], "cancer": row["cancer"], "platform": row["platform"]}
    return adata


def standardize_result(res: pd.DataFrame, row: Dict[str, str], lr_pathway: Dict[Tuple[str, str], str]) -> List[Dict[str, object]]:
    out = []
    if res is None or res.empty:
        return out
    for _, r in res.iterrows():
        ligand = str(r.get("ligand", "")).upper()
        receptor = str(r.get("receptor", "")).upper()
        source = str(r.get("source", ""))
        target = str(r.get("target", ""))
        pathway = lr_pathway.get((ligand, receptor), "")
        out.append(
            {
                "method": METHOD,
                "dataset": row["dataset"],
                "sample_id": row["sample_id"],
                "cancer": row["cancer"],
                "platform": row["platform"],
                "sender_compartment": source,
                "receiver_compartment": target,
                "ligand": ligand,
                "receptor": receptor,
                "pathway": pathway,
                "lr_means": r.get("lr_means", ""),
                "cellphone_pvals": r.get("cellphone_pvals", ""),
                "ligand_props": r.get("ligand_props", ""),
                "receptor_props": r.get("receptor_props", ""),
            }
        )
    return out


def run_sample(row: Dict[str, str], annot_row: Dict[str, str], lr_df: pd.DataFrame, args: argparse.Namespace):
    adata = build_adata(row, annot_row, lr_df, args.min_total_counts, args.min_cells)
    available = set(adata.var_names)
    sample_lr = lr_df[lr_df["ligand"].isin(available) & lr_df["receptor"].isin(available)].copy().reset_index(drop=True)
    if sample_lr.empty:
        raise ValueError(f"No pilot LR pairs available in expression matrix for {row['sample_id']}")
    resource = sample_lr[["ligand", "receptor"]].copy()
    li.method.cellphonedb(
        adata,
        groupby="spot_compartment",
        resource=resource,
        expr_prop=args.expr_prop,
        min_cells=args.min_cells,
        n_perms=args.n_perms,
        seed=args.seed,
        n_jobs=args.n_jobs,
        use_raw=False,
        key_added="liana_cellphonedb_res",
        verbose=False,
    )
    res = adata.uns.get("liana_cellphonedb_res", pd.DataFrame())
    lr_pathway = {(r["ligand"], r["receptor"]): r["pathway"] for _, r in sample_lr.iterrows()}
    candidates = standardize_result(res, row, lr_pathway)
    summary = {
        "method": METHOD,
        "dataset": row["dataset"],
        "sample_id": row["sample_id"],
        "cancer": row["cancer"],
        "platform": row["platform"],
        "spots_used": int(adata.n_obs),
        "genes_used": int(adata.n_vars),
        "compartments_used": int(adata.obs["spot_compartment"].nunique()),
        "lr_pairs_requested": int(lr_df.shape[0]),
        "lr_pairs_available": int(sample_lr.shape[0]),
        "candidate_rows": int(len(candidates)),
        "expr_prop": args.expr_prop,
        "min_cells": args.min_cells,
        "n_perms": args.n_perms,
    }
    return candidates, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run LIANA CellPhoneDB-style pilot baseline.")
    parser.add_argument("--max-samples", type=int, default=1)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--output-prefix", default="liana_cellphonedb_pilot")
    parser.add_argument("--expr-prop", type=float, default=0.05)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--min-total-counts", type=int, default=1)
    parser.add_argument("--n-perms", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--n-jobs", type=int, default=2)
    args = parser.parse_args()

    manifest_rows = read_table(PROJECT / "docs/minimal_input_manifest.tsv")
    annot_rows = {row["sample_id"]: row for row in read_table(PROJECT / "docs/spot_annotation_manifest.tsv")}
    lr_df = read_lr_panel(PROJECT / "docs/pilot_lr_panel.tsv")
    selected = choose_samples(manifest_rows, args.sample_id, args.max_samples)
    out_dir = PROJECT / "results/task_c_liana_baseline"
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
        except Exception as exc:
            failures.append({"dataset": row.get("dataset", ""), "sample_id": row.get("sample_id", ""), "error_type": type(exc).__name__, "error_message": str(exc)})

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
        "outputs": {"candidates": rel(cand_path), "sample_summary": rel(summ_path), "failures": rel(fail_path)},
        "parameters": {"max_samples": args.max_samples, "sample_id": args.sample_id, "output_prefix": args.output_prefix, "expr_prop": args.expr_prop, "min_cells": args.min_cells, "n_perms": args.n_perms, "seed": args.seed, "n_jobs": args.n_jobs},
        "note": "LIANA CellPhoneDB-style pilot baseline using coarse computational spot compartments; not a final benchmark result.",
    }
    json_path.write_text(json.dumps(summary_json, indent=2, sort_keys=True) + "\n")
    manifest_path = run_dir / "outputs_manifest.tsv"
    with manifest_path.open("w") as handle:
        handle.write("path\ttype\tdescription\n")
        handle.write(f"{rel(cand_path)}\ttsv\tStandardized LIANA CellPhoneDB pilot candidates\n")
        handle.write(f"{rel(summ_path)}\ttsv\tSample-level LIANA pilot summary\n")
        handle.write(f"{rel(fail_path)}\ttsv\tLIANA pilot failures, if any\n")
        handle.write(f"{rel(json_path)}\tjson\tLIANA pilot run summary\n")
        handle.write(f"{rel(manifest_path)}\ttsv\tRun output manifest\n")
    print(json.dumps(summary_json, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
