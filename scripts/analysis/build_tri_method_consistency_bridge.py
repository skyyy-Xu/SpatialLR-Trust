#!/usr/bin/env python3
"""Build tri-method consistency support across stdlib, COMMOT, and LIANA pilot baselines."""
from __future__ import annotations
import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "manual")
KEY_COLS = ["dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment", "ligand", "receptor", "pathway"]
METHODS = ["stdlib_marker_lr_pilot", "commot_pilot_v1", "liana_cellphonedb_pilot"]

def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT))

def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: List[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def key(row: Dict[str, str]) -> Tuple[str, ...]:
    values = []
    for col in KEY_COLS:
        val = row.get(col, "")
        if col in {"ligand", "receptor"}:
            val = val.upper()
        values.append(val)
    return tuple(values)

def float_or_zero(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0

def index_rows(rows: List[Dict[str, str]]) -> Dict[Tuple[str, ...], Dict[str, str]]:
    out: Dict[Tuple[str, ...], Dict[str, str]] = {}
    for row in rows:
        out.setdefault(key(row), row)
    return out

def support_label(row: Dict[str, object]) -> str:
    parts = []
    if row["stdlib_supported"]:
        parts.append("stdlib")
    if row["commot_supported"]:
        parts.append("commot")
    if row["liana_supported"]:
        parts.append("liana")
    return "+".join(parts)

def bucket_label(row: Dict[str, object]) -> str:
    label = support_label(row).replace("+", "_")
    return label + ("_only" if row["method_support_count"] == 1 else "")

def summarize(counter: Counter) -> Dict[str, object]:
    union = int(counter["union"])
    at_least_two = int(counter["at_least_two"])
    all_three = int(counter["stdlib_commot_liana"])
    return {
        "union": union,
        "all_three": all_three,
        "at_least_two": at_least_two,
        "stdlib_commot_only": int(counter["stdlib_commot"]),
        "stdlib_liana_only": int(counter["stdlib_liana"]),
        "commot_liana_only": int(counter["commot_liana"]),
        "stdlib_only": int(counter["stdlib_only"]),
        "commot_only": int(counter["commot_only"]),
        "liana_only": int(counter["liana_only"]),
        "stdlib": int(counter["stdlib"]),
        "commot": int(counter["commot"]),
        "liana": int(counter["liana"]),
        "all_three_fraction_of_union": all_three / union if union else 0.0,
        "at_least_two_fraction_of_union": at_least_two / union if union else 0.0,
    }

def main() -> None:
    stdlib_path = PROJECT / "results/task_c_pilot_baseline/stdlib_marker_lr_pilot_candidates.tsv"
    commot_path = PROJECT / "results/task_c_commot_baseline/commot_all_candidates.tsv"
    liana_path = PROJECT / "results/task_c_liana_baseline/liana_cellphonedb_all_candidates.tsv"
    out_dir = PROJECT / "results/task_c_method_consistency"
    run_dir = PROJECT / "runs" / RUN_ID
    out_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    stdlib_rows = read_tsv(stdlib_path)
    commot_rows = read_tsv(commot_path)
    liana_rows = read_tsv(liana_path)
    stdlib = index_rows(stdlib_rows)
    commot = index_rows(commot_rows)
    liana = index_rows(liana_rows)
    pair_rows = []
    for k in sorted(set(stdlib) | set(commot) | set(liana)):
        std = stdlib.get(k)
        com = commot.get(k)
        lia = liana.get(k)
        base = dict(zip(KEY_COLS, k))
        std_s = 1 if std else 0
        com_s = 1 if com else 0
        lia_s = 1 if lia else 0
        ext_count = com_s + lia_s
        method_count = std_s + com_s + lia_s
        row = {
            **base,
            "method_support_count": method_count,
            "method_support_label": "",
            "stdlib_supported": std_s,
            "commot_supported": com_s,
            "liana_supported": lia_s,
            "external_method_support_count": ext_count,
            "external_method_support_fraction": ext_count / 2.0,
            "stdlib_pilot_score": float_or_zero(std.get("pilot_score", "0")) if std else 0.0,
            "stdlib_adjacency_enrichment": float_or_zero(std.get("adjacency_enrichment", "0")) if std else 0.0,
            "commot_mean_score": float_or_zero(com.get("commot_mean_score", "0")) if com else 0.0,
            "liana_lr_means": float_or_zero(lia.get("lr_means", "0")) if lia else 0.0,
            "liana_cellphone_pvals": float_or_zero(lia.get("cellphone_pvals", "1")) if lia else 1.0,
            "sender_spots": (com or std or {}).get("sender_spots", ""),
            "receiver_spots": (com or std or {}).get("receiver_spots", ""),
        }
        row["method_support_label"] = support_label(row)
        pair_rows.append(row)
    by_dataset = defaultdict(Counter)
    by_pathway = defaultdict(Counter)
    by_sample = defaultdict(Counter)
    for row in pair_rows:
        label = bucket_label(row)
        buckets = (by_dataset[row["dataset"]], by_pathway[row["pathway"]], by_sample[(row["dataset"], row["sample_id"], row["cancer"])])
        for bucket in buckets:
            bucket["union"] += 1
            bucket[label] += 1
            bucket["stdlib"] += int(row["stdlib_supported"])
            bucket["commot"] += int(row["commot_supported"])
            bucket["liana"] += int(row["liana_supported"])
            if int(row["method_support_count"]) >= 2:
                bucket["at_least_two"] += 1
    pair_path = out_dir / "tri_method_consistency_pair_support.tsv"
    dataset_path = out_dir / "tri_method_consistency_by_dataset.tsv"
    pathway_path = out_dir / "tri_method_consistency_by_pathway.tsv"
    sample_path = out_dir / "tri_method_consistency_by_sample.tsv"
    summary_path = out_dir / "tri_method_consistency_summary.json"
    manifest_path = run_dir / "outputs_manifest.tsv"
    pair_fields = KEY_COLS + ["method_support_count", "method_support_label", "stdlib_supported", "commot_supported", "liana_supported", "external_method_support_count", "external_method_support_fraction", "stdlib_pilot_score", "stdlib_adjacency_enrichment", "commot_mean_score", "liana_lr_means", "liana_cellphone_pvals", "sender_spots", "receiver_spots"]
    summary_fields = ["union", "all_three", "at_least_two", "stdlib_commot_only", "stdlib_liana_only", "commot_liana_only", "stdlib_only", "commot_only", "liana_only", "stdlib", "commot", "liana", "all_three_fraction_of_union", "at_least_two_fraction_of_union"]
    write_tsv(pair_path, pair_rows, pair_fields)
    write_tsv(dataset_path, [{"dataset": d, **summarize(c)} for d, c in sorted(by_dataset.items())], ["dataset"] + summary_fields)
    write_tsv(pathway_path, [{"pathway": p, **summarize(c)} for p, c in sorted(by_pathway.items())], ["pathway"] + summary_fields)
    write_tsv(sample_path, [{"dataset": d, "sample_id": s, "cancer": ca, **summarize(c)} for (d, s, ca), c in sorted(by_sample.items())], ["dataset", "sample_id", "cancer"] + summary_fields)
    total = Counter()
    for row in pair_rows:
        total["union"] += 1
        total[bucket_label(row)] += 1
        total["stdlib"] += int(row["stdlib_supported"])
        total["commot"] += int(row["commot_supported"])
        total["liana"] += int(row["liana_supported"])
        if int(row["method_support_count"]) >= 2:
            total["at_least_two"] += 1
    summary = {
        "run_id": RUN_ID,
        "methods": METHODS,
        "stdlib_rows": len(stdlib_rows),
        "commot_rows": len(commot_rows),
        "liana_rows": len(liana_rows),
        **summarize(total),
        "outputs": {"pair_support": rel(pair_path), "by_dataset": rel(dataset_path), "by_pathway": rel(pathway_path), "by_sample": rel(sample_path)},
        "note": "Tri-method consistency bridge across stdlib pilot, COMMOT pilot-parameter baseline, and LIANA CellPhoneDB-style baseline; not final method benchmark evidence.",
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    with manifest_path.open("w") as handle:
        handle.write("path\ttype\tdescription\n")
        for path, typ, desc in [(pair_path, "tsv", "Per-candidate method support across stdlib, COMMOT, and LIANA"), (dataset_path, "tsv", "Tri-method support counts by dataset"), (pathway_path, "tsv", "Tri-method support counts by pathway"), (sample_path, "tsv", "Tri-method support counts by sample"), (summary_path, "json", "Tri-method consistency summary"), (manifest_path, "tsv", "Run output manifest")]:
            handle.write(f"{rel(path)}\t{typ}\t{desc}\n")
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
