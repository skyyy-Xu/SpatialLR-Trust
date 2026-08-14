# Task D Fake LR Null Summary

Date: 2026-07-09

## Result

Run `20260709_1511_fake-lr-null-pilot` attached an expression-matched fake ligand-receptor pair null to the `stdlib_marker_lr_pilot` candidate table.

- Candidate rows: 4225
- Fake draws requested per candidate: 100
- Candidates without usable fake pairs: 9
- Empirical p <= 0.05: 1239
- Empirical p <= 0.10: 1531

## Outputs

- Null table: `results/task_d_null_models/stdlib_marker_lr_fake_pair_null.tsv`
- Summary JSON: `results/task_d_null_models/stdlib_marker_lr_fake_pair_null_summary.json`

## Interpretation

This null tests whether an observed candidate score is stronger than expression-nearest fake ligand/receptor combinations drawn from the same pilot LR gene pool. It is still attached to the pilot baseline and does not replace full external communication-tool inference.
