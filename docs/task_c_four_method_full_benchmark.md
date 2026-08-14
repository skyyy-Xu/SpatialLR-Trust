# Four-Method Full-Sample Benchmark

The full benchmark compares the stdlib marker baseline, COMMOT, LIANA and full-expression CellChat across all 34 MVP
samples. LR-component support is keyed by dataset, sample, sender, receiver, ligand and receptor without pathway, so
database-specific pathway labels do not create artificial pair mismatches. Pathway-direction support is reported as a
separate output.

The LR-component union contains 108,187 keys: 101,761 supported by one method, 2,836 by two, 2,793 by three and 797 by
all four. CellChat contributes 102,306 unique LR-component keys from 93,790 standardized rows. Most single-method keys
reflect the deliberate scope mismatch between full `CellChatDB.human` and the curated pilot panel used by the other
three methods, rather than a calibrated false-positive rate.

Restricting LR identities to the pilot-comparable panel yields 7,346 keys: 920 supported by one method, 2,836 by two,
2,793 by three and 797 by all four; 6,426 receive at least two-method support. At the pathway-direction level, the union
contains 25,305 keys, including 4,579 supported by at least two methods and 690 by all four.

These outputs strengthen the method-consistency evidence layer but do not by themselves change the v3 trust score.
A CellChat-aware score should be evaluated as a separate sensitivity analysis because the candidate databases and
pathway naming systems are not identical.
