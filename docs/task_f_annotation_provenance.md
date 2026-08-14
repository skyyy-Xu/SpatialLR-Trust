# Task F GSE280634 Annotation Provenance

## Scope

F2 retains two annotation layers for the GSE280634 Xenium object:

- `native_label`: the deposited author-provided `Cell Type` value.
- `harmonized_compartment`: an explicit coarse mapping used only for cross-platform comparison.

The harmonized label is an analytical abstraction, not a replacement for the native annotation and not an experimental validation of cell identity.

## Mapping policy

The complete mapping is frozen in `docs/task_f_gse280634_annotation_mapping.tsv`. Conversion must stop if a native label is absent from that table. No substring, regular-expression or nearest-name fallback is allowed.

The target compartments are `tumor/epithelial`, `T/NK`, `B/plasma`, `myeloid`, `CAF/fibroblast`, `endothelial` and `other/unknown`. Smooth muscle is conservatively retained as `other/unknown` because the deposited label does not establish a fibroblast or CAF identity. The author residual label `Other` also remains `other/unknown`.

## Expression and coordinates

- Counts: deposited `layers/counts`, promoted to canonical `X` without a duplicate count layer.
- Coordinates: deposited `obsm/spatial`, retained as canonical `obsm/spatial`.
- Sample identity: deposited `batch`, retained as `batch_native` and namespaced to `GSE280634_batch_<value>` in `sample_id`.
- Cell identity: deposited unique observation index, retained as the canonical observation index.

## Filtering and QC

No arbitrary transcript-count cutoff is introduced in F2. Cells with a zero sum in `layers/counts` are excluded and reported. Duplicate cell IDs, non-finite coordinates, missing sample/native labels or unmapped labels are hard failures rather than silent filters. All 280 targeted genes are retained unless the conversion fails because feature identifiers are not unique.

F2 establishes a reproducible input contract. It does not show that downstream ligand-receptor calls are true, causal or experimentally validated.
