# CellChat v4 Score Sensitivity

## Decision

No formal SpatialLR-Trust v4 score is published from this analysis. The current formal pilot score remains v3. The most defensible immediate use of CellChat is a descriptive `CellChat-confirmed v3 high` subset containing 159 of the 218 v3 high candidates.

## Purpose

The completed full-expression CellChat benchmark provides a third external communication method, but adding a method to a composite score can change both weights and tier gates. This analysis quantifies those changes before any score-version release.

The input is the unchanged 4,225-row v3 candidate universe with CellChat diagnostics. No output is written under `scores/`.

## Scenarios

Eight scenarios were evaluated:

1. The unchanged v3 reference.
2. A confirmation-only rule that preserves v3 scores but requires COMMOT, LIANA and CellChat for high confidence.
3. Six reweighted grids combining consistency weights of 0.10, 0.20 or 0.30 with either an all-three or at-least-two external-method high gate.

For reweighted scenarios, the non-method evidence terms retain their relative proportions and are rescaled to fill `1 - consistency_weight`. COMMOT, LIANA and CellChat contribute equally within the consistency term. The high and medium score thresholds remain 0.75 and 0.55, and every high tier must pass all three null gates.

## Results

| Scenario | High | Medium | Low | Tier changes | High retained from v3 | High with CellChat | Spearman vs v3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| v3 reference | 218 | 1,976 | 2,031 | 0 | 100.00% | 72.94% | 1.000 |
| CellChat confirmation only | 159 | 2,035 | 2,031 | 59 | 72.94% | 100.00% | 1.000 |
| Weight 0.10, all three | 158 | 1,642 | 2,425 | 480 | 72.48% | 100.00% | 0.980 |
| Weight 0.10, at least two | 215 | 1,585 | 2,425 | 423 | 98.62% | 73.49% | 0.980 |
| Weight 0.20, all three | 159 | 1,674 | 2,392 | 424 | 72.94% | 100.00% | 0.991 |
| Weight 0.20, at least two | 216 | 1,617 | 2,392 | 367 | 99.08% | 73.61% | 0.991 |
| Weight 0.30, all three | 159 | 1,760 | 2,306 | 462 | 72.94% | 100.00% | 0.984 |
| Weight 0.30, at least two | 216 | 1,703 | 2,306 | 405 | 99.08% | 73.61% | 0.984 |

The confirmation-only scenario is the smallest perturbation: it demotes 59 v3 high candidates to medium and leaves every other tier unchanged. The equal-three 0.20 scenario has the best rank stability among reweighted grids (Spearman 0.991), but changes 424 tiers: 59 high to medium, 363 medium to low and two low to medium.

At-least-two gates preserve 215–216 high candidates but only 73–74% of those high candidates have CellChat support. These gates therefore mostly reproduce the v3 high set rather than using CellChat as independent confirmation.

## Why v4 is deferred

- CellChat and the other methods share ligand-receptor database assumptions and are not independent ground truth.
- CellChat used its full human database, whereas COMMOT and LIANA were run on the curated pilot LR panel; equal support weights do not remove this coverage asymmetry.
- VEGF was a clear CellChat failure case in the preceding diagnostic (0/291 v3 candidates supported), so a strict all-three gate may remove biologically plausible candidates for database or method-specific reasons.
- The grid is not calibrated against downstream target activation or an independent positive-reference set.

Before a formal v4 release, the candidate formula and gate should be evaluated against pathway target evidence, curated public positive references and the documented VEGF failure case. Until then, v3 remains the score, and the 159-candidate CellChat-confirmed high subset is descriptive sensitivity evidence only.

## Outputs

- `results/task_e_cellchat_v4_sensitivity/cellchat_v4_sensitivity_scenarios.tsv`
- `results/task_e_cellchat_v4_sensitivity/cellchat_v4_sensitivity_transitions.tsv`
- `results/task_e_cellchat_v4_sensitivity/cellchat_v4_sensitivity_candidates.tsv.gz`
- `results/task_e_cellchat_v4_sensitivity/cellchat_v4_sensitivity_summary.json`

