# Software Environment Summary
Date: 2026-08-08
Run IDs: `20260802_1040_runtime-environment-lock`; `20260808_1342_runtime-reconstruction-v8`
Status: lightweight environment and command-entry summary for the current pilot benchmark.

## Scope
This document indexes the project-local Python and R environments, key package versions, the versioned runtime lock and manuscript-facing command entry points. It does not reinstall software or rerun long analyses.

## Runtime
| item | value | note |
|---|---|---|
| project_root | $PROJECT_ROOT | Project-root placeholder |
| platform | Linux-6.12.0-211.16.1.el10_2.0.1.x86_64-x86_64-with-glibc2.39 | Detected on the xulab login node |
| default_python | Python 3.12.13 | Used by stdlib-only scripts |
| project_python | Python 3.12.13 | envs/spatiallr-py312-commot/bin/python |
| project_r | Rscript (R) version 4.3.3 (2024-02-29) | envs/spatiallr-r-cellchat-cytosignal/bin/Rscript |
| git | git version 2.52.0 | Version-control tool |
| runtime_lock | configs/runtime_environment_lock_v1.json | Versioned installed-runtime evidence and reconstruction boundaries |
| runtime_reconstruction_v8 | results/task_f_reproducibility/runtime_reconstruction/20260808_1342_runtime-reconstruction-v8/release_validation.json | Same-cluster empty-environment release evidence |

## Key Python Packages
| import | distribution | version | role |
|---|---|---|---|
| commot | commot | 0.0.3 | COMMOT spatial communication baseline |
| liana | liana | 1.8.0 | LIANA CellPhoneDB-style baseline |
| scanpy | scanpy | 1.11.5 | AnnData/scRNA-seq and spatial data utilities |
| anndata | anndata | 0.12.19 | Annotated matrix container |
| numpy | numpy | 1.26.4 | Numerical arrays |
| scipy | scipy | 1.17.1 | Scientific computing |
| pandas | pandas | 2.3.3 | Tabular data processing |
| matplotlib | matplotlib | 3.10.9 | Figure generation |
| seaborn | seaborn | 0.13.2 | Statistical plotting |
| networkx | networkx | 3.6.1 | Graph utilities |
| ot | pot | 0.9.7 | Optimal transport backend used by COMMOT |
| sklearn | scikit-learn | 1.9.0 | Machine-learning utilities |
| pyarrow | pyarrow | 19.0.1 | Columnar file I/O used by post-F8 audits |
| rdata | rdata | 1.1.0 | Read-only RData inspection support |
| shapely | shapely | 2.1.2 | Geometry validation and synthetic policy |
| xarray | xarray | 2026.4.0 | Array dependency used by rdata |

## Key R Packages
| package | version | role |
|---|---|---|
| CellChat | 2.2.0.9001 | Full-expression CellChat baseline |
| NMF | 0.28 | Pinned local-source CellChat dependency |
| Biobase | 2.62.0 | Bioconductor expression container support |
| ComplexHeatmap | 2.18.0 | CellChat visualization dependency |
| ggplot2 | 3.5.2 | R plotting |
| igraph | 2.1.4 | CellChat graph utilities |

## Command Entry Points
| command id | command | role |
|---|---|---|
| minimal_input_qc | `python3 scripts/server/full_minimal_input_qc.py` | Minimal input QC |
| annotation | `python3 scripts/server/marker_compartment_annotation.py --mode all --run-id <RUN_ID>` | Coarse marker compartment annotation |
| marker_baseline | `python3 scripts/server/stdlib_marker_lr_pilot.py` | Marker-baseline LR candidate generation |
| commot_baseline | `$PY scripts/server/run_commot_pilot_baseline.py --max-samples 999 --distance-threshold 2.5 --cot-nitermax 2000` | COMMOT baseline |
| liana_baseline | `RUN_ID=<RUN_ID> MAX_SAMPLES=999 OUTPUT_PREFIX=liana_cellphonedb_all sbatch scripts/server/run_liana_cellphonedb_baseline.sh` | LIANA CellPhoneDB-style baseline |
| cellchat_baseline | `RUN_ID=<RUN_ID> sbatch scripts/server/run_cellchat_full_sample_baseline.sh` | Full-expression CellChat baseline |
| cellchat_v3_diagnostics | `python3 scripts/server/build_cellchat_v3_consistency_diagnostics.py` | Post-freeze CellChat-v3 consistency |
| spatial_null | `python3 scripts/server/spatial_null_for_pilot_lr.py --permutations 50 --seed 20260709` | Spatial-coordinate null |
| label_null | `python3 scripts/server/label_null_for_pilot_lr.py --permutations 30 --seed 20260709` | Label-permutation null |
| fake_lr_null | `python3 scripts/server/fake_lr_null_for_pilot.py --draws 100 --seed 20260709` | Expression-matched fake LR null |
| v3_score | `python3 scripts/server/build_pilot_trust_score_v3.py` | Pilot v3 score |
| target_proxy | `python3 scripts/server/build_target_activation_layer.py --max-samples 999 --output-prefix target_activation_all` | Receiver target-expression proxy |
| public_evidence | `TOP_N=40 python3 scripts/server/build_public_evidence_check_top40.py` | Public evidence categorization |
| recurrence_summary | `python3 scripts/server/build_recurrence_by_axis_summary.py --run-id <RUN_ID>` | Axis-level recurrence summary |
| repro_checks | `RUN_ID=<RUN_ID> scripts/server/run_repro_checks.sh` | Path and numeric/text invariant checks |
| runtime_lock_check | `python3 scripts/server/check_runtime_environment_lock.py` | Installed Python/R runtime lock validation |
| runtime_reconstruction_v8_release | `python3 scripts/server/check_runtime_reconstruction_v8_release.py` | Independent bounded reconstruction release gate |
| bilingual_check | `python3 scripts/server/check_bilingual_manuscript.py` | English-Chinese protected-content parity |

## Boundary
This summary is an index, not the lock itself. `docs/runtime_environment_lock.md`, `configs/runtime_environment_lock_v1.json` and the v8 release evidence distinguish installed snapshots, passed same-cluster empty-environment reconstruction, unresolved artifact-byte identity, and unvalidated container or cross-platform claims.

## Outputs
- `results/task_f_reproducibility/software_environment.tsv`
- `results/task_f_reproducibility/python_package_versions.tsv`
- `results/task_f_reproducibility/r_package_versions.tsv`
- `results/task_f_reproducibility/command_entrypoints.tsv`
- `results/task_f_reproducibility/software_environment_summary.json`
