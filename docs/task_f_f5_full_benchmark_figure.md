# Task F F5 Full Benchmark Figure

Status: `PASS_F5_FULL_BENCHMARK_FIGURE`

## Figure contract

Core conclusion: the frozen score has a modest common-universe discrimination
gain and lower non-low FDP, but unresolved hard negatives and weaker
technical-replicate stability limit the strength of the claim.

Archetype: quantitative grid with a hero comparison panel.

Evidence chain:

- panel a: common-universe AUROC and average precision with section-cluster
  bootstrap intervals
- panel b: false-positive rates for four hard-negative families
- panel c: complete-series monotonic response to spike and corruption
- panel d: technical-replicate rank stability across the three sections

Python matplotlib is the exclusive plotting and export backend. Planned
exports are editable SVG, PDF, 600-dpi TIFF and PNG preview under
`figures/task_f_f5_full_benchmark`.

## Statistical and integrity notes

- biological replicates: three GSE243168 sections
- technical replicates: five simulation replicates per section
- bootstrap: 2,000 draws clustered by section, seed 20260726
- hard-negative FPR denominator: 60 rows per family
- monotonicity series: 60 per scorer except LIANA, which has 30 because of
  frozen resource scope
- rank stability: ten pairwise replicate comparisons per section and scorer
- no hypothesis-test p-values or multiple-comparison claims
- all quantitative panels trace to
  `results/task_f_cross_platform/f5_semisim_full_score`

## Validated export

Final figure job `47802` passed programmatic and visual QA:

- editable-text SVG
- PDF
- 600-dpi TIFF at 4,946 x 3,990 pixels
- lossless LZW TIFF size: approximately 2.0 MB
- PNG preview
- zero export-validation failures
- no overlapping labels or panels in the final PNG inspection

The first job failed before publication while serializing Pillow DPI metadata.
The next export passed programmatic QA, but manual inspection found panel c
labels too close to panel d. The final export increased lower-panel spacing and
shortened panel c labels, then losslessly compressed the TIFF without changing
data, pixels, resolution or statistics. The final publication pass also
removed Matplotlib SVG line-end whitespace without changing vector semantics.

The figure intentionally retains the mixed result: panel a shows the modest
common-universe point-estimate gain, while panels b-d expose unresolved hard
negatives, coarse monotonic response and weaker technical-replicate stability.
