# Pilot SpatialLR-Trust Scoring Schema v3
Date: 2026-07-09
Run ID: `20260709_2040_pilot-score-v3`

## Scope
This schema extends v2 by adding LIANA CellPhoneDB-style support as a second external-method consistency feature.
The scored universe remains the 4,225 stdlib pilot candidates because current null-model evidence was generated for
that candidate table.

## Inputs
- Base score: `scores/spatiallr_trust_score_pilot_v2.tsv`
- Tri-method support: `results/task_c_method_consistency/tri_method_consistency_pair_support.tsv`

## Formula
```text
external_method_support_fraction = (COMMOT same-key support + LIANA same-key support) / 2
spatiallr_trust_score_pilot_v3 = 0.40 * null_support_mean
                               + 0.15 * pilot_score_percentile
                               + 0.15 * recurrence_support
                               + 0.10 * annotation_quality
                               + 0.20 * external_method_support_fraction
```

## Tier Rules
- `high_pilot_v3`: score >= 0.75, all three pilot null p-values <= 0.10, and both COMMOT and LIANA support the same key.
- `medium_pilot_v3`: score >= 0.55 but not high.
- `low_pilot_v3`: score < 0.55.
Empty null p-values are treated conservatively as unsupported for high-tier gating.

## Results
- Rows scored: 4,225
- COMMOT-supported rows: 3,301
- LIANA-supported rows: 3,448
- Both external methods supported: 2,924
- One external method supported: 901
- No external method supported: 400
- `high_pilot_v3`: 218
- `medium_pilot_v3`: 1,976
- `low_pilot_v3`: 2,031

For comparison, v2 tiers were 222 high, 1,986 medium, and 2,017 low. The v3 high tier is slightly stricter because it
requires support from both external baselines, not just COMMOT.

## Outputs
- `scores/spatiallr_trust_score_pilot_v3.tsv`
- `results/task_e_scores/spatiallr_trust_score_pilot_v3_summary.json`

## Interpretation Boundary
This is still a pilot score. It improves cross-method consistency evidence, but final scoring should add stronger
cell-type annotation, spatially explicit external baselines where possible, pathway target evidence, and benchmark-level
positive/negative references.
