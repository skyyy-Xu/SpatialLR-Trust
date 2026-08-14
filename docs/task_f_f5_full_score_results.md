# Task F F5 Full Score Results

Status: `PASS_F5_SEMISIM_FULL_SCORE_MIXED_PERFORMANCE`

## Validated result

The frozen TF15 formula was applied without fitting or retuning to all 1,140
F5 truth-axis rows. Independent validation passed:

- 1,080 resource rows with scores
- 60 fixed-fake rows retained as unscored out-of-scope controls
- 76 condition-axis recurrence keys
- three biological sections
- five technical simulation replicates per section
- eight primary/scoped metric rows
- 16 section-cluster bootstrap summaries
- zero validation failures

The score SHA-256 is
`cf5dfd8039a0bb6a3ede2032b4417436d48f5ced59cfbc88baa04c907be2bf60`.

Tier counts were 81 high, 554 medium, 445 low and 60 fixed-fake
out-of-scope. These remain computational score strata rather than calibrated
truth probabilities.

## Discrimination

On the 510-row common-three-method universe (390 positives and 120 hard
negatives), SpatialLR-Trust had the highest point estimates:

| Scorer | AUROC | Section-bootstrap 95% interval | Average precision |
| --- | ---: | ---: | ---: |
| SpatialLR-Trust | 0.621 | 0.526-0.666 | 0.824 |
| COMMOT | 0.583 | 0.579-0.589 | 0.809 |
| LIANA | 0.553 | 0.498-0.608 | 0.813 |
| CellChat | 0.494 | 0.463-0.494 | 0.755 |

The bootstrap unit is the biological section and all five technical
replicates remain clustered. With only three sections, interval comparisons
are descriptive and do not establish statistical superiority.

On the 1,020-row all-resource labelled universe, SpatialLR-Trust had AUROC
0.549 and average precision 0.781. COMMOT covers the same labelled universe
and had AUROC 0.550 and average precision 0.769. The score therefore does not
show a uniform discrimination gain across universes.

## False discoveries

Empirical false-discovery proportions use only positive/negative labelled
resource rows. Fixed fake and unmodified-reference rows are excluded:

| Selection | Selected | True positive | False positive | Empirical FDP |
| --- | ---: | ---: | ---: | ---: |
| Trust non-low | 613 | 496 | 117 | 0.191 |
| COMMOT detected | 1,001 | 779 | 222 | 0.222 |
| LIANA detected | 505 | 390 | 115 | 0.228 |
| CellChat detected | 912 | 702 | 210 | 0.230 |
| Trust high | 81 | 57 | 24 | 0.296 |

The non-low score selection reduces FDP relative to raw method detections, but
the high tier does not. High-tier false positives are concentrated in
reporter-absent (14/60) and spatially separated (10/60) hard negatives. The
high tier must therefore not be described as a validated high-confidence
class.

## Robustness limits

SpatialLR-Trust was fully monotonic in only 31.7% of spike-strength series,
compared with 100% for COMMOT and LIANA and 93.3% for CellChat. Its
full-series monotonic fractions were 33.3% for coordinate corruption and
31.7-33.3% for label and expression corruption.

Mean pairwise technical-replicate Spearman correlations for SpatialLR-Trust
were 0.814, 0.788 and 0.643 across the three sections. Raw-method values were
generally 0.974-0.999. The discrete null-pass and recurrence components improve
some negative filtering but introduce ranking instability and coarse
responses.

## Interpretation

The retained result is mixed:

- common-universe discrimination and non-low empirical FDP improve
  descriptively;
- all-resource AUROC does not improve over COMMOT;
- the high tier fails important hard negatives;
- monotonic response and technical-replicate rank stability are weaker than
  the raw methods.

No weight, threshold, recurrence rule, tier or metric universe is changed
after inspection. The result supports a benchmark narrative about the
trade-off between conservative evidence integration and score granularity,
not a claim of universal method superiority.

The next analysis should focus on the reporter-absent and spatial-separation
failure mechanisms and on whether continuous null evidence can be evaluated
as a predeclared future score variant in a new development/validation split.
It must not overwrite the frozen F5 result.
