# SpatialLR-Trust

SpatialLR-Trust is an auditable **operational prioritization and failure-analysis workflow** for spatial ligand-receptor candidates. It harmonizes candidate calls from multiple inference routes and evaluates marker-baseline candidates with finite spatial-coordinate, label-permutation and expression-matched fake-pair screens, descriptive recurrence, annotation-design priors and external-method consistency.

> **Release status:** this is a public review snapshot derived from source commit `413563b6c2e86abbac6286d2a8c3b2fdd2db31e8`. It is not yet a versioned software release. Author metadata, the final licence, the external archive DOI and independent off-cluster execution remain pending corresponding-author approval.

## Current evidence

- Public pilot benchmark: 34 tumour spatial-transcriptomics samples across three GEO datasets.
- Frozen marker-baseline universe: 4,225 scored candidates; 218 in the high-priority pilot tier.
- Negative controls: spatial-coordinate permutation, label permutation and expression-matched fake ligand-receptor pairs.
- Baselines: marker route, COMMOT, LIANA/CellPhoneDB-style calls and a post-freeze CellChat sensitivity analysis.
- Stage 4 sensitivity: leave-one-biological-sample-out recurrence changed the pilot high tier from 218 to 217; candidate-gene exclusion changed 178 target-proxy support states and reduced high-tier support from 134 to 128.
- Truth-labelled semi-simulation: mixed discrimination and robustness; the workflow is not claimed to be universally superior or probabilistically calibrated.

Operational tiers quantify robustness to the specified computational screens. They do **not** establish pathway activation, resolved sender-receiver directionality, experimental signalling truth or cross-platform performance.

## Repository contents

```text
configs/                 compact runtime metadata
data/                    source-data availability notes (no raw GEO files)
docs/                    method, data and reproduction documentation
figures/                 final PNG/PDF/SVG figures (no TIFF files)
results/                 selected compact derived benchmark outputs
scores/                  frozen pilot v3 score table
scripts/analysis/        core analysis and figure scripts
scripts/release/         deterministic release-snapshot validator
```

## Verify this snapshot

Only Python 3.10+ from the standard library is needed for the release integrity gate:

```bash
python scripts/release/check_release.py
```

The gate checks the release manifest, forbidden/private paths, file-size policy, Python syntax and frozen benchmark invariants. Full regeneration requires the public source datasets and method-specific Python/R environments described in `docs/REPRODUCIBILITY.md`.

## Data

Raw GEO data are not redistributed. Dataset accessions, source links and analytic roles are listed in `docs/data_inventory.tsv` and `docs/DATA_AVAILABILITY.md`. Selected compact derived tables are included to make manuscript claims inspectable.

## Citation and licence

Citation metadata and licensing remain proposals for author review. See `CITATION.cff.template`, `LICENSE_PROPOSAL.md` and `docs/RELEASE_REVIEW_CHECKLIST.md`. The repository is publicly viewable, but no reuse licence is granted until the owners approve and add the final licence files.

## Source

- Public repository: https://github.com/skyyy-Xu/SpatialLR-Trust
- Authoritative internal source commit: `413563b6c2e86abbac6286d2a8c3b2fdd2db31e8`
- Snapshot date: 2026-08-15
