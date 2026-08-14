# Runtime Environment Lock

Run IDs: `20260802_1040_runtime-environment-lock`; `20260808_1342_runtime-reconstruction-v8`

## Scope

This lock records the installed Python and R runtimes used by the project. It is an auditable installed snapshot, not a claim of cross-platform or bitwise reconstruction.

## Locked evidence

- Python Python 3.12.13 with 125 pinned installed distributions and a passing `pip check`.
- R 4.3.3 on `x86_64-conda-linux-gnu` with 259 installed packages.
- A Linux-64 explicit conda base lock for the R environment.
- CellChat and NMF local source archives with independently checked sizes and SHA-256 values.
- A versioned JSON manifest and per-artifact hash table.

## Reconstruction boundary

- Python is an installed package snapshot without a complete distribution-hash lock.
- The R conda base has an explicit Linux-64 lock, while CellChat and NMF are source-installed and locked separately.
- Same-cluster empty-environment reconstruction passed under one committed source lock: the Python branch reproduced the exact 125-distribution installed-version set, and the R branch reproduced the exact Linux-64 URL base plus six canonical source overlays and the 259-package final inventory.
- This reconstruction does not establish Python wheel or sdist byte identity, conda package-artifact byte identity, cross-platform reconstruction or a validated container.
- No Docker or Apptainer image has been built or validated.
- CytoSignal remains unresolved and is not installed or claimed.

## Validation

```bash
python3 scripts/server/check_runtime_environment_lock.py
```

The installed-runtime validator is read-only with respect to installed environments. It checks snapshot hashes, current package state, source archive locks and forbidden absolute-path or credential markers. The v8 reconstruction evidence additionally passed preflight (14/14), Python (34/34), R (55/55), compute aggregate (34/34) and independent login-node verification (10/10) gates. The last is the workflow's internal release gate, not public-release validation.
