# CellChat-v3 Consistency Diagnostics

## Purpose

This analysis asks whether full-expression CellChat support is enriched among candidates that already received higher SpatialLR-Trust v3 confidence. It is a descriptive sensitivity analysis. It does not alter the v3 score formula, score values or confidence tiers.

## Join definition

- Input score table: `scores/spatiallr_trust_score_pilot_v3.tsv` (4,225 candidates).
- CellChat evidence: `results/task_c_four_method_full_benchmark/four_method_lr_component_support.tsv.gz`.
- Matching key: dataset, sample, cancer, sender compartment, receiver compartment, normalized ligand and normalized receptor.
- Pathway is excluded from the matching key because database pathway labels are not directly comparable across methods.
- All 4,225 v3 candidates matched the four-method support table; no candidate was dropped.

## Main result

CellChat supported 815 of 4,225 v3 candidates (19.29%). Support was strongly ordered by the pre-existing v3 tier:

| v3 tier | Candidates | CellChat supported | Support fraction |
| --- | ---: | ---: | ---: |
| High | 218 | 159 | 72.94% |
| Medium | 1,976 | 552 | 27.94% |
| Low | 2,031 | 104 | 5.12% |

The CellChat support fraction was 14.24-fold higher in the high tier than in the low tier. This ordering is consistent with the v3 score capturing cross-method reliability despite being defined before full-expression CellChat was added.

## Null-model consistency

Among 225 candidates passing all three null gates (`spatial_null_p`, `label_null_p` and `fake_lr_null_p` at 0.10), 159 were CellChat-supported (70.67%). Among the remaining 4,000 candidates, 656 were CellChat-supported (16.40%). Thus, CellChat support was 4.31-fold enriched among all-null-pass candidates.

All 218 high-tier candidates passed the three null gates by construction; 159 were additionally supported by CellChat. Seven medium-tier candidates passed all null gates, but none had CellChat support. No low-tier candidate passed all three null gates.

## Dataset and pathway patterns

CellChat support varied across dataset-cancer groups: melanoma 38.40% (134/349), colorectal cancer 21.05% (225/1,069), NSCLC 17.10% (445/2,602) and DLBCL 5.37% (11/205). This variation should not be read as a biological ranking because sample composition, compartment annotation and database coverage differ.

The largest supported pathway groups were TGFB (174), SPP1 (170), MIF (100) and NOTCH (97). MIF, TGFB and SPP1 had support fractions of 31.06%, 30.37% and 28.91%, respectively. VEGF was a clear failure case: none of 291 v3 VEGF candidates matched CellChat support under the pathway-agnostic LR-component key. This may reflect database/LR identity coverage, compartment granularity or method sensitivity and should be investigated before any v4 scoring change.

## Interpretation boundary

The monotonic tier pattern is computational agreement, not experimental validation. CellChat is not an independent ground truth and shares ligand-receptor database assumptions with other communication tools. These diagnostics therefore support retaining the v3 ranking for now and motivate a later, separately calibrated v4 sensitivity analysis; they do not justify adding CellChat directly to the formal score without threshold and failure-case evaluation.

## Outputs

- `results/task_e_cellchat_v3_diagnostics/cellchat_v3_candidate_diagnostics.tsv.gz`
- `results/task_e_cellchat_v3_diagnostics/cellchat_support_by_v3_tier.tsv`
- `results/task_e_cellchat_v3_diagnostics/cellchat_support_by_pathway.tsv`
- `results/task_e_cellchat_v3_diagnostics/cellchat_support_by_dataset.tsv`
- `results/task_e_cellchat_v3_diagnostics/cellchat_support_by_null_gate.tsv`
- `results/task_e_cellchat_v3_diagnostics/cellchat_v3_diagnostics_summary.json`

