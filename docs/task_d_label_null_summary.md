# Task D Label Null Summary

Date: 2026-07-09

## Result

Run `20260709_1137_label-null-pilot` attached a spot-compartment label permutation null to the `stdlib_marker_lr_pilot` candidate table.

- Candidate rows: 4225
- Samples: 34
- Permutations: 30
- Empirical p <= 0.05: 915
- Empirical p <= 0.10: 1070

| Dataset | Candidate rows |
| --- | ---: |
| GSE283052 | 1069 |
| GSE292299 | 2807 |
| GSE300445 | 349 |

## Outputs

- Null table: `results/task_d_null_models/stdlib_marker_lr_label_null.tsv`
- Summary JSON: `results/task_d_null_models/stdlib_marker_lr_label_null_summary.json`

## Interpretation

This null tests whether the observed candidate score depends on the observed coarse compartment labels rather than a random redistribution of labels across fixed tissue coordinates. It is still tied to the pilot marker-compartment baseline and does not replace external communication-tool statistics.
