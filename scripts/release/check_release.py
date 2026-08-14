#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import py_compile
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "release_manifest.tsv"
FORBIDDEN_SUFFIXES = {".tiff", ".docx", ".h5", ".h5ad", ".rds", ".tar", ".zip"}
FORBIDDEN_TEXT = [
    "/home/" + "user/xiaotian",
    "/Users/" + "skyyy",
    "BEGIN OPENSSH " + "PRIVATE KEY",
]
MAX_BYTES = 50 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, label: str, checks: list[dict]) -> None:
    checks.append({"check": label, "status": "PASS" if condition else "FAIL"})
    if not condition:
        raise AssertionError(label)


def main() -> int:
    checks = []
    required = [
        ROOT / "README.md",
        ROOT / "scores/spatiallr_trust_score_pilot_v3.tsv",
        ROOT / "scripts/analysis/build_pilot_trust_score_v3.py",
        ROOT / "figures/task_e_global_results_overview/spatiallr_global_results_overview.png",
    ]
    require(all(path.is_file() for path in required), "required release files", checks)

    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) >= 100, "manifest has at least 100 curated files", checks)
    for row in rows:
        path = ROOT / row["path"]
        require(path.is_file(), f"manifest path exists: {row['path']}", checks)
        require(path.stat().st_size == int(row["bytes"]), f"size matches: {row['path']}", checks)
        require(sha256(path) == row["sha256"], f"hash matches: {row['path']}", checks)

    files = [path for path in ROOT.rglob("*") if path.is_file() and ".git" not in path.parts]
    require(all(path.stat().st_size < MAX_BYTES for path in files), "no file reaches 50 MiB", checks)
    require(all(path.suffix.lower() not in FORBIDDEN_SUFFIXES for path in files), "forbidden binary types absent", checks)

    text_failures = []
    for path in files:
        if path.suffix not in {".py", ".R", ".sh", ".md", ".tsv", ".json", ".txt", ".yml", ".yaml", ".cff"}:
            continue
        text = path.read_text(encoding="utf-8")
        if any(pattern in text for pattern in FORBIDDEN_TEXT):
            text_failures.append(str(path.relative_to(ROOT)))
    require(not text_failures, "private paths and key material absent", checks)

    score_path = ROOT / "scores/spatiallr_trust_score_pilot_v3.tsv"
    with score_path.open(encoding="utf-8", newline="") as handle:
        score_rows = list(csv.DictReader(handle, delimiter="\t"))
    tiers = Counter(row["confidence_tier_pilot_v3"] for row in score_rows)
    require(len(score_rows) == 4225, "pilot score has 4,225 rows", checks)
    require(tiers == {"high_pilot_v3": 218, "medium_pilot_v3": 1976, "low_pilot_v3": 2031}, "pilot tier counts match", checks)

    summary_path = ROOT / (
        "results/task_f_reproducibility/tf48_stage4_revision_analyses/"
        "20260810_1631_tf48-stage4-revision-analyses/summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require(summary["pilot_original_high"] == 218 and summary["pilot_loo_high"] == 217, "pilot LOO summary matches", checks)
    require(summary["target"]["support_flips"] == 178, "candidate-exclusion flips match", checks)
    require(summary["target"]["high_excluded_supported_rows"] == 128, "candidate-excluded high support matches", checks)
    require(summary["f5_original_high"] == 81 and summary["f5_loo_high"] == 41, "semi-simulation LOO summary matches", checks)

    python_files = sorted((ROOT / "scripts/analysis").glob("*.py")) + [Path(__file__)]
    for path in python_files:
        py_compile.compile(str(path), doraise=True)
    require(len(python_files) >= 15, "core Python scripts compile", checks)

    report = {
        "status": "PASS_SPATIALLR_GITHUB_REVIEW_SNAPSHOT",
        "checks": len(checks),
        "passed": sum(item["status"] == "PASS" for item in checks),
        "failed": sum(item["status"] == "FAIL" for item in checks),
        "details": checks,
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "FAIL_SPATIALLR_GITHUB_REVIEW_SNAPSHOT", "error": str(exc)}, indent=2), file=sys.stderr)
        raise
