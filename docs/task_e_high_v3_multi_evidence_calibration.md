# High-v3 Multi-evidence Calibration

## Purpose

This analysis cross-calibrates three evidence layers within the unchanged 218 SpatialLR-Trust v3 high candidates:

- full-expression CellChat same-key support;
- receiver-compartment target-expression support;
- the existing manually reviewed public-evidence top40 matrix.

The goal is to define a computational validation shortlist and identify discordant failure cases. It does not modify v3 or publish a v4 score.

## Main result

| CellChat | Target proxy | High-v3 rows | Interpretation |
| --- | --- | ---: | --- |
| Supported | Supported | 108 | Joint computational consensus |
| Supported | Not supported | 51 | CellChat-only; inspect target-panel sensitivity |
| Not supported | Supported | 26 | Target-only; inspect CellChat/database coverage |
| Not supported | Not supported | 33 | Highest-priority unresolved or filtered cases |

The 108 joint rows represent 49.54% of all v3 high candidates, 67.92% of CellChat-supported high candidates and 80.60% of target-supported high candidates.

CellChat and target support were positively associated within the high tier (odds ratio 2.69; two-sided Fisher exact `p=0.00170`). This agreement is stronger than expected from the marginal support rates, but it does not establish biological truth because both layers are computational and may share expression-related biases.

## Dataset stability

Joint support was similar across the four dataset-cancer groups:

| Dataset | Cancer | High v3 | Joint support | Fraction |
| --- | --- | ---: | ---: | ---: |
| GSE283052 | Colorectal cancer | 75 | 39 | 52.00% |
| GSE292299 | DLBCL | 13 | 7 | 53.85% |
| GSE292299 | NSCLC | 101 | 47 | 46.53% |
| GSE300445 | Melanoma | 29 | 15 | 51.72% |

The narrow 46.5–53.8% range argues against the joint subset being driven by one cohort alone, although all cohorts still use coarse compartment annotations.

## Pathway patterns and failure cases

ECM-integrin had the highest joint fraction (18/23, 78.26%), followed by EGF/ERBB (6/9, 66.67%) and MIF (26/41, 63.41%). SPP1 contributed the largest absolute joint count (28), followed by MIF (26) and ECM-integrin (18).

VEGF remained the clearest failure case: all five v3 high VEGF candidates lacked both CellChat and target-proxy support. This is not enough to label them false positives because the CellChat database and the curated target panel can both miss valid VEGF biology. These candidates should be audited against LR identity coverage, receiver target-panel coverage and original dataset evidence before score rules are changed.

The 77 discordant candidates (51 CellChat-only and 26 target-only) are retained explicitly. They are method-diagnostic cases rather than rows to discard automatically.

## Public evidence overlay

All 40 candidates in the existing public-evidence review matched the high-v3 matrix. Thirty-eight were jointly supported by CellChat and the target proxy, while two were target-only. Evidence grades were:

- 7 direct-axis, same-cancer rows;
- 21 direct-axis, other-cancer rows;
- 8 direct ECM-receptor, other-cancer rows;
- 4 broad SPP1-integrin rows.

This concentration is expected because the top40 was selected from target-supported v3 high candidates. It is useful for prioritization but cannot be treated as an unbiased validation rate. The original-paper audit found zero same-sample candidate-specific validation rows.

## Decision boundary

The 108 joint rows form the strongest current computational shortlist. They may be described as `v3 high with CellChat and target-expression support`. They must not be called experimentally validated interactions.

Formal v4 scoring remains deferred. The next evidence work should focus on the 77 discordant candidates and five VEGF neither-supported high candidates, and should seek independent pathway activation or original-paper validation rather than adding another uncalibrated score term.

## Outputs

- `results/task_e_high_v3_multi_evidence/high_v3_multi_evidence_matrix.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_joint_cellchat_target_consensus.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_cellchat_target_discordant.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_multi_evidence_by_dataset.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_multi_evidence_by_pathway.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_multi_evidence_categories.tsv`
- `results/task_e_high_v3_multi_evidence/high_v3_multi_evidence_summary.json`

