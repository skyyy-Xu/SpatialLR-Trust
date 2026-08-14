# Reproducibility

## Level 1: release integrity

Run from the repository root:

```bash
python scripts/release/check_release.py
```

This checks file hashes, file policy, private-path exclusion, Python syntax and the frozen numeric summary. It does not rerun biological analyses.

## Level 2: selected output regeneration

Set the project root explicitly if invoking a script from another directory:

```bash
export PROJECT="$PWD"
```

The core order is:

1. marker-based compartment annotation;
2. marker ligand-receptor candidate generation;
3. COMMOT, LIANA and CellChat baseline calls;
4. candidate-key harmonization;
5. spatial, label and fake-pair null screens;
6. frozen v3 operational score;
7. target-proxy and recurrence sensitivity analyses;
8. semi-simulation benchmark and figures.

Scripts are provided under `scripts/analysis/`. They preserve the validated project logic but use a portable project-root default. Raw GEO matrices are intentionally absent.

## Level 3: full benchmark regeneration

Full regeneration requires downloading the public source datasets in `docs/data_inventory.tsv`, reconstructing method-specific Python and R environments, and recording each long computation through a scheduler or equivalent batch system. The internal same-cluster reconstruction gate passed, but no Docker/Apptainer image or independent off-cluster clean execution has yet been validated.

## Frozen interpretation boundary

The included score is an operational prioritization score, not a calibrated probability. The semi-simulation contains only three biological sections with clustered technical replicates. Cross-platform routes that failed input-eligibility gates did not produce performance evidence.
