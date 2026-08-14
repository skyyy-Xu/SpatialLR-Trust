#!/usr/bin/env python3
"""Build all-sample four-method support with full-expression CellChat results."""
from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
RUN_ID = "20260713_1445_four-method-full-benchmark"
METHODS = ("stdlib", "commot", "liana", "cellchat")
LR_COLS = ("dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment", "ligand", "receptor")
PATHWAY_COLS = LR_COLS[:5] + ("pathway",)


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def norm(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.strip().upper())


def compartment(value: str) -> str:
    return value.strip().lower()


def receptor_components(row: dict[str, str]) -> list[str]:
    ligand = norm(row.get("ligand", ""))
    tokens = [norm(token) for token in row.get("interaction_name", "").split("_") if norm(token)]
    if len(tokens) >= 2 and tokens[0] == ligand:
        return list(dict.fromkeys(tokens[1:]))
    receptor = norm(row.get("receptor", ""))
    return [receptor] if receptor else []


def adapt_existing(rows: list[dict[str, str]], method: str) -> list[dict[str, str]]:
    return [{**row, "_method": method} for row in rows]


def adapt_cellchat(rows: list[dict[str, str]], sample_cancer: dict[str, str]) -> list[dict[str, str]]:
    adapted = []
    for row in rows:
        for receptor in receptor_components(row):
            adapted.append(
                {
                    **row,
                    "cancer": sample_cancer[row["sample_id"]],
                    "sender_compartment": row["sender"],
                    "receiver_compartment": row["receiver"],
                    "receptor": receptor,
                    "cellchat_original_receptor": row.get("receptor", ""),
                    "_method": "cellchat",
                }
            )
    return adapted


def lr_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["dataset"], row["sample_id"], row["cancer"].lower(), compartment(row["sender_compartment"]),
        compartment(row["receiver_compartment"]), norm(row["ligand"]), norm(row["receptor"]),
    )


def pathway_key(row: dict[str, str]) -> tuple[str, ...]:
    key = lr_key(row)
    return key[:5] + (norm(row.get("pathway", "")),)


def index(rows: list[dict[str, str]], key_fn) -> dict[tuple[str, ...], dict[str, str]]:
    result = {}
    for row in rows:
        result.setdefault(key_fn(row), row)
    return result


def support_label(flags: dict[str, int]) -> str:
    return "+".join(method for method in METHODS if flags[method])


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    distribution = Counter(int(row["method_support_count"]) for row in rows)
    comparable = [row for row in rows if int(row["pilot_lr_pair"]) == 1]
    comparable_distribution = Counter(int(row["method_support_count"]) for row in comparable)
    return {
        "union": len(rows),
        "support_count_distribution": {str(key): distribution[key] for key in sorted(distribution)},
        "at_least_two": sum(count for key, count in distribution.items() if key >= 2),
        "all_four": distribution[4],
        "pilot_comparable_union": len(comparable),
        "pilot_comparable_support_count_distribution": {
            str(key): comparable_distribution[key] for key in sorted(comparable_distribution)
        },
        "pilot_comparable_at_least_two": sum(count for key, count in comparable_distribution.items() if key >= 2),
        "pilot_comparable_all_four": comparable_distribution[4],
    }


def main() -> None:
    raw = {
        "stdlib": read_tsv(PROJECT / "results/task_c_pilot_baseline/stdlib_marker_lr_pilot_candidates.tsv"),
        "commot": read_tsv(PROJECT / "results/task_c_commot_baseline/commot_all_candidates.tsv"),
        "liana": read_tsv(PROJECT / "results/task_c_liana_baseline/liana_cellphonedb_all_candidates.tsv"),
        "cellchat": read_tsv(PROJECT / "results/task_c_cellchat_all_sample_baseline/cellchat_all_candidates.tsv.gz"),
    }
    sample_manifest = read_tsv(
        PROJECT / "results/task_c_cellchat_all_sample_baseline/cellchat_full_sample_manifest.tsv"
    )
    sample_cancer = {row["sample_id"]: row["cancer"] for row in sample_manifest}
    adapted = {
        "stdlib": adapt_existing(raw["stdlib"], "stdlib"),
        "commot": adapt_existing(raw["commot"], "commot"),
        "liana": adapt_existing(raw["liana"], "liana"),
        "cellchat": adapt_cellchat(raw["cellchat"], sample_cancer),
    }
    lr_index = {method: index(adapted[method], lr_key) for method in METHODS}
    pilot_pairs = {
        (norm(row["ligand"]), norm(row["receptor"]))
        for method in ("stdlib", "commot", "liana")
        for row in adapted[method]
    }
    all_keys = sorted(set().union(*(set(lr_index[method]) for method in METHODS)))
    pair_rows: list[dict[str, object]] = []
    for key in all_keys:
        flags = {method: int(key in lr_index[method]) for method in METHODS}
        method_pathways = {
            method: norm(lr_index[method].get(key, {}).get("pathway", "")) for method in METHODS
        }
        cellchat = lr_index["cellchat"].get(key, {})
        pair_rows.append(
            {
                **dict(zip(LR_COLS, key)),
                "pilot_lr_pair": int((key[5], key[6]) in pilot_pairs),
                "method_support_count": sum(flags.values()),
                "method_support_label": support_label(flags),
                **{f"{method}_supported": flags[method] for method in METHODS},
                **{f"{method}_pathway": method_pathways[method] for method in METHODS},
                "cellchat_probability": cellchat.get("probability", ""),
                "cellchat_p_value": cellchat.get("p_value", ""),
                "cellchat_interaction_name": cellchat.get("interaction_name", ""),
                "cellchat_original_receptor": cellchat.get("cellchat_original_receptor", ""),
            }
        )
    pathway_index = {method: index(adapted[method], pathway_key) for method in METHODS}
    pathway_keys = sorted(set().union(*(set(pathway_index[method]) for method in METHODS)))
    pathway_rows: list[dict[str, object]] = []
    for key in pathway_keys:
        flags = {method: int(key in pathway_index[method]) for method in METHODS}
        pathway_rows.append(
            {
                **dict(zip(PATHWAY_COLS, key)),
                "method_support_count": sum(flags.values()),
                "method_support_label": support_label(flags),
                **{f"{method}_supported": flags[method] for method in METHODS},
            }
        )
    out_dir = PROJECT / "results/task_c_four_method_full_benchmark"
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_path = out_dir / "four_method_lr_component_support.tsv.gz"
    pathway_path = out_dir / "four_method_pathway_direction_support.tsv.gz"
    sample_path = out_dir / "four_method_support_by_sample.tsv"
    dataset_path = out_dir / "four_method_support_by_dataset.tsv"
    pair_fields = list(LR_COLS) + [
        "pilot_lr_pair", "method_support_count", "method_support_label", "stdlib_supported", "commot_supported",
        "liana_supported", "cellchat_supported", "stdlib_pathway", "commot_pathway", "liana_pathway",
        "cellchat_pathway", "cellchat_probability", "cellchat_p_value", "cellchat_interaction_name",
        "cellchat_original_receptor",
    ]
    pathway_fields = list(PATHWAY_COLS) + [
        "method_support_count", "method_support_label", "stdlib_supported", "commot_supported", "liana_supported",
        "cellchat_supported",
    ]
    write_tsv(pair_path, pair_rows, pair_fields)
    write_tsv(pathway_path, pathway_rows, pathway_fields)
    by_sample: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    by_dataset: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        by_sample[(str(row["dataset"]), str(row["sample_id"]), str(row["cancer"]))].append(row)
        by_dataset[str(row["dataset"])].append(row)
    summary_fields = [
        "union", "at_least_two", "all_four", "pilot_comparable_union", "pilot_comparable_at_least_two",
        "pilot_comparable_all_four",
    ]
    sample_rows = [
        {"dataset": key[0], "sample_id": key[1], "cancer": key[2], **{name: summarize(rows)[name] for name in summary_fields}}
        for key, rows in sorted(by_sample.items())
    ]
    dataset_rows = [
        {"dataset": key, **{name: summarize(rows)[name] for name in summary_fields}}
        for key, rows in sorted(by_dataset.items())
    ]
    write_tsv(sample_path, sample_rows, ["dataset", "sample_id", "cancer"] + summary_fields)
    write_tsv(dataset_path, dataset_rows, ["dataset"] + summary_fields)
    overall = summarize(pair_rows)
    pathway_distribution = Counter(int(row["method_support_count"]) for row in pathway_rows)
    summary = {
        "run_id": RUN_ID,
        "methods": list(METHODS),
        "raw_rows": {method: len(raw[method]) for method in METHODS},
        "unique_lr_component_keys": {method: len(lr_index[method]) for method in METHODS},
        **overall,
        "pathway_direction_union": len(pathway_rows),
        "pathway_direction_support_count_distribution": {
            str(key): pathway_distribution[key] for key in sorted(pathway_distribution)
        },
        "pathway_direction_at_least_two": sum(count for key, count in pathway_distribution.items() if key >= 2),
        "pathway_direction_all_four": pathway_distribution[4],
        "samples": len(sample_rows),
        "datasets": len(dataset_rows),
        "key_definition": "LR-component support excludes pathway from the key; pathway-direction support is reported separately.",
        "comparability_boundary": "Stdlib, COMMOT and LIANA used the curated pilot LR panel, while CellChat used all CellChatDB.human categories. pilot_comparable metrics restrict LR identities to the pilot panel.",
        "outputs": {
            "lr_component_support": str(pair_path.relative_to(PROJECT)),
            "pathway_direction_support": str(pathway_path.relative_to(PROJECT)),
            "by_sample": str(sample_path.relative_to(PROJECT)),
            "by_dataset": str(dataset_path.relative_to(PROJECT)),
        },
    }
    (out_dir / "four_method_full_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
