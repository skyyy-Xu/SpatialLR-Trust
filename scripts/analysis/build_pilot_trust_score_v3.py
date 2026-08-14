#!/usr/bin/env python3
"""Build v3 pilot SpatialLR-Trust score with COMMOT and LIANA support."""
from __future__ import annotations
import csv
import json
import os
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT = Path(os.environ.get("PROJECT", str(Path(__file__).resolve().parents[2])))
RUN_ID = os.environ.get("RUN_ID", "manual")
BASE_SCORE = PROJECT / "scores/spatiallr_trust_score_pilot_v2.tsv"
CONSISTENCY = PROJECT / "results/task_c_method_consistency/tri_method_consistency_pair_support.tsv"
OUT = PROJECT / "scores/spatiallr_trust_score_pilot_v3.tsv"
SUMMARY = PROJECT / "results/task_e_scores/spatiallr_trust_score_pilot_v3_summary.json"
RUN_DIR = PROJECT / "runs" / RUN_ID
KEY_FIELDS = ["dataset", "sample_id", "cancer", "sender_compartment", "receiver_compartment", "ligand", "receptor", "pathway"]

def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def key(row: Dict[str, str]) -> Tuple[str, ...]:
    out = []
    for field in KEY_FIELDS:
        val = row.get(field, "")
        if field in {"ligand", "receptor"}:
            val = val.upper()
        out.append(val)
    return tuple(out)

def f(row: Dict[str, str], field: str, default: float = 0.0) -> float:
    try:
        return float(row.get(field, default))
    except Exception:
        return default

def tier(row: Dict[str, object]) -> str:
    score = float(row["spatiallr_trust_score_pilot_v3"])
    null_ok = all(f(row, p, 1.0) <= 0.10 for p in ["spatial_null_p", "label_null_p", "fake_lr_null_p"])
    both_external = int(row["external_method_support_count"]) == 2
    if score >= 0.75 and null_ok and both_external:
        return "high_pilot_v3"
    if score >= 0.55:
        return "medium_pilot_v3"
    return "low_pilot_v3"

def main() -> None:
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    base_rows = read_tsv(BASE_SCORE)
    support_rows = read_tsv(CONSISTENCY)
    support = {key(row): row for row in support_rows}
    out_rows = []
    for row in base_rows:
        sup = support.get(key(row), {})
        commot = int(float(sup.get("commot_supported", row.get("commot_supported", 0)) or 0))
        liana = int(float(sup.get("liana_supported", 0) or 0))
        external_count = commot + liana
        external_fraction = external_count / 2.0
        score = 0.40 * f(row, "null_support_mean") + 0.15 * f(row, "pilot_score_percentile") + 0.15 * f(row, "recurrence_support") + 0.10 * f(row, "annotation_quality") + 0.20 * external_fraction
        new = dict(row)
        new.update({
            "tri_method_support_count": int(float(sup.get("method_support_count", 1) or 1)),
            "commot_supported_v3": commot,
            "liana_supported": liana,
            "external_method_support_count": external_count,
            "external_method_support_fraction": external_fraction,
            "liana_lr_means": f(sup, "liana_lr_means"),
            "liana_cellphone_pvals": f(sup, "liana_cellphone_pvals", 1.0),
            "spatiallr_trust_score_pilot_v3": round(score, 6),
            "notes_v3": "v3 adds LIANA CellPhoneDB-style support as a second external-method consistency feature; null evidence remains stdlib-pilot based.",
        })
        new["confidence_tier_pilot_v3"] = tier(new)
        out_rows.append(new)
    fieldnames = list(base_rows[0].keys()) + ["tri_method_support_count", "commot_supported_v3", "liana_supported", "external_method_support_count", "external_method_support_fraction", "liana_lr_means", "liana_cellphone_pvals", "spatiallr_trust_score_pilot_v3", "confidence_tier_pilot_v3", "notes_v3"]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)
    tiers = Counter(row["confidence_tier_pilot_v3"] for row in out_rows)
    summary = {
        "run_id": RUN_ID,
        "rows_scored": len(out_rows),
        "commot_supported_rows": sum(int(row["commot_supported_v3"]) for row in out_rows),
        "liana_supported_rows": sum(int(row["liana_supported"]) for row in out_rows),
        "both_external_supported_rows": sum(1 for row in out_rows if int(row["external_method_support_count"]) == 2),
        "one_external_supported_rows": sum(1 for row in out_rows if int(row["external_method_support_count"]) == 1),
        "no_external_supported_rows": sum(1 for row in out_rows if int(row["external_method_support_count"]) == 0),
        "tiers": dict(sorted(tiers.items())),
        "formula": "0.40*null_support_mean + 0.15*pilot_score_percentile + 0.15*recurrence_support + 0.10*annotation_quality + 0.20*((COMMOT+LIANA)/2)",
        "high_tier_rule": "score>=0.75, all three pilot null p-values<=0.10, and both COMMOT and LIANA support the same key",
        "outputs": {"score_table": str(OUT.relative_to(PROJECT)), "summary": str(SUMMARY.relative_to(PROJECT))},
        "note": "Pilot v3 score remains scoped to stdlib-scored candidates because null models were generated for that universe.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest = RUN_DIR / "outputs_manifest.tsv"
    with manifest.open("w") as handle:
        handle.write("path\ttype\tdescription\n")
        handle.write(f"{OUT.relative_to(PROJECT)}\ttsv\tPilot SpatialLR-Trust v3 score with COMMOT and LIANA support\n")
        handle.write(f"{SUMMARY.relative_to(PROJECT)}\tjson\tPilot v3 score summary\n")
        handle.write(f"{manifest.relative_to(PROJECT)}\ttsv\tRun output manifest\n")
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
