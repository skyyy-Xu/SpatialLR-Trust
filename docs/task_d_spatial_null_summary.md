# Task D Spatial Null Summary

Date: 2026-07-09

## Result

Run `20260709_0115_spatial-null-pilot` attached a spatial coordinate permutation null to the `stdlib_marker_lr_pilot` candidate table.

- Candidate rows: 4225
- Samples: 34
- Permutations: 50
- Empirical p <= 0.05: 1711
- Empirical p <= 0.10: 1783

| Dataset | Candidate rows |
| --- | ---: |
| GSE283052 | 1069 |
| GSE292299 | 2807 |
| GSE300445 | 349 |

## Outputs

- Null table: `results/task_d_null_models/stdlib_marker_lr_spatial_null.tsv`
- Summary JSON: `results/task_d_null_models/stdlib_marker_lr_spatial_null_summary.json`

## Interpretation

This is the first candidate-level empirical null interface. It evaluates whether observed compartment adjacency support exceeds randomized spatial coordinate assignments under the pilot marker-compartment labels. It does not replace external communication tools or full statistical inference.
