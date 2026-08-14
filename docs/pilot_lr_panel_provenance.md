# Pilot Ligand-Receptor Panel Provenance

Date: 2026-08-02
Run ID: `20260802_1120_manuscript-presubmission-review`
Status: retrospective provenance and claim-boundary record for the frozen pilot panel.

## Origin

`docs/pilot_lr_panel.tsv` was added in commit
`af88001ed18dfeddef8bc61ed8ba19c28f9070db` as part of the initial marker
ligand-receptor pilot baseline, before the benchmark candidate outputs, null
models and v3 tiers were generated. The file contains 15 directed
ligand-receptor pairs spanning TGFB, chemokine, VEGF, SPP1, MIF, NOTCH, WNT,
EGF/ERBB and ECM-integrin pathways.

## Intended use

The panel was defined for pipeline development, pathway diversity and interface
testing across the marker baseline, COMMOT and LIANA routes. It fixes the formal
v3 candidate universe on which all three pilot null models were calculated.

## Boundary

- The panel is project-defined; it is not claimed to be a complete or
  systematically curated ligand-receptor resource.
- It is not an independently validated positive set or biological ground truth.
- Reference 6 documents signalling-resource context but does not define these
  exact 15 pairs.
- Full-expression CellChat used the broader `CellChatDB.human` resource and is
  therefore retained as post-freeze descriptive evidence rather than a formal
  v3 input.
- Adding or removing pairs would define a new candidate universe and would
  require new null models and a new score version; this review does not alter
  `docs/pilot_lr_panel.tsv` or any frozen result.
