# Running This Pipeline

This document gives the exact order to run every script in this repository,
from a fresh clone to the final failure-mode review. It exists specifically
so a reviewer or collaborator does not have to reconstruct the pipeline order
from `PROJECT_HANDOVER.md`'s narrative.

**Read this first:** this pipeline was built and run entirely on a single
16GB-RAM Linux workstation and takes multiple days of wall-clock time in
total (see per-stage estimates below). It is not designed to run end-to-end
in one sitting. Every stage below is checkpointed / resumable by design —
rerunning a stage's script after a partial failure will skip already-completed
work, not redo it from scratch.

---

## 0. Environment setup

```bash
# R 4.3.3 required. Install renv if not already present:
install.packages("renv")

# Restore the exact pinned package versions:
renv::restore()
```

**SymSim requires a manual, non-renv patch** (a `BiocGenerics::rank()` S4
compatibility conflict that requires source-level patching — runtime patching
does not work because SymSim is bytecode-compiled):

```bash
Rscript code/01_simulation/install_symsim_patched.R
```

This clones `YosefLab/SymSim` at its current HEAD (not pinned to a specific
commit — a known, documented reproducibility gap; see `PATCH_NOTES.md`).

**Storage note:** the full pipeline produces roughly 340GB of intermediate
and processed data (raw simulated matrices, per-method embeddings). This
repository tracks only the small, essential artifacts (ground truth, final
result tables, code); the large intermediate data is expected to live on
local/external storage, not in git. See `README.md`'s Data section.

---

## 1. Simulation generation and calibration (Steps 1.1–1.8)

Run in `code/01_simulation/`, in this order:

1. `calibrate_splatter.R`, `calibrate_splatter_bcv.R`, `calibrate_splatter_bcv2.R`, `calibrate_splatter_depth_dropout.R` — Splatter parameter calibration.
2. `calibrate_scdesign3.R` — PBMC 3k reference prep and scDesign3 calibration.
3. `param_dict.R`, then `param_grid.R` — builds the invariant 10,940-row parameter grid (`data/simulated/param_grid.csv`).
4. `simulate_splatter.R`, `simulate_scdesign3.R`, `simulate_symsim.R` — main-grid generation, one simulator at a time. **Do not use** the files suffixed `_v1_INVALID_sparsity_inert.R` — these are preserved buggy originals, kept as documented evidence, not meant to be run.
5. `simulate_null_control_splatter.R`, `simulate_null_control_scdesign3.R`, `simulate_null_control_symsim.R` — null-control set (45 files total).
6. `validate_output_inventory.R` — full validation pass against calibration tables.
7. `extract_ground_truth_splatter.R`, `extract_ground_truth_scdesign3.R`, `extract_ground_truth_symsim.R` — ground-truth signal-subspace extraction (outputs tracked in `data/simulated/ground_truth/`).
8. `convert_to_sce.R` — standardizes everything into `SingleCellExperiment` format.

**Estimated runtime:** multiple days in aggregate (SymSim's larger fits alone can take 1.5–2+ hours each; ~10,940 files × 3 simulators).

**Known caveat:** 6 specific Splatter `run_id`s (3284, 3285, 6565, 6569, 9849, 9850) are permanently excluded after exhaustive root-cause investigation found no resolvable cause — this is expected and documented, not a bug to chase.

---

## 2. Real-data curation (Step 2)

Run in `code/02_real_data/`:

1. `download_pbmc68k.R`, `download_pancreas_datasets.R` (Baron, Muraro, Segerstolpe), `download_tabula_sapiens_lung.R`, `download_tasic2018_brain.R`.
2. `qc_filter_one_dataset.R` (via `run_qc_all_datasets.sh`).
3. `harmonize_one_dataset.R` (via `run_harmonize_all_datasets.sh`).
4. `reannotate_pbmc68k.R` — corrects PBMC68k's original labeling gap.
5. `verify_tabula_sapiens_lung.R`, `verify_tasic2018_brain.R`.
6. `build_dataset_inventory.R` — produces `data/real_data_inventory.csv`.

**Requires internet access** to the original data repositories (GEO, `scRNAseq`/`TENxPBMCData` Bioconductor packages).

---

## 3. EDA Checkpoint 1 (simulated-output verification)

**No standalone script exists for this checkpoint** — it was run as ad hoc,
interactive R commands during a working session and never saved as a
reusable file. See `docs/eda_checkpoint1_reconstruction_note.md` for the
exact logic to reconstruct it, and its tracked outputs in
`results/step1_eda/`.

---

## 4. EDA Checkpoint 2 (real-dataset technical characterization)

Run in `code/eda_checkpoint2/`, in this order:

`technical_summary_one_dataset.R` → `pancreas_batch_pca.R` →
`pancreas_batch_pca_log1p.R` → `pancreas_batch_pca_per_dataset.R` →
`compute_simulated_phase_space.R` → `plot_phase_space_overlay.R` →
`regime_alignment_analysis.R`.

**Known caveat:** `irlba()` calls in this checkpoint were run without a
preceding `set.seed()` — a documented reproducibility gap, not yet fixed.

---

## 5. Preprocessing and dimensionality reduction (Step 3)

Run in `code/03_preprocessing/`, one method at a time (each is independent
and can run in any order relative to the others, but each must complete
fully before Step 4 for that method):

1. `step3_1_raw_pca_full.R` (Raw PCA)
2. `step3_2_libnorm_pca_full.R` (Library-normalized PCA)
3. `step3_3_log_pca_full.R` (Log-PCA)
4. `step3_4_shiftedlog_pca_full.R` (Shifted-Log PCA)
5. `step3_5_sctransform_with_loadings.R` (Pearson-residual/SCTransform v2 — **use this final version**, not `step3_5_sctransform_pca_full.R`, which is the superseded original)
6. `step3_6_glmpca_real.R` + `step3_6_glmpca_simulated.R` (GLM-PCA)
7. `step3_loadings_backfill.R` — backfills analytically-reconstructed gene loadings for methods 1–4.
8. `step3_10_manifest.R` — builds the master `embedding_manifest.csv` integrity scan and manifest (**required by every downstream Step 4–6 script**).

**Estimated runtime:** several hours per method, GLM-PCA and SCTransform v2 being the slowest.

---

## 6. EDA Checkpoint 3 (representative embedding visualization)

```bash
Rscript code/eda_checkpoint3/eda_checkpoint3.R
```
Requires Step 3.10's `embedding_manifest.csv` to exist first.

---

## 7. Subspace-recovery metrics (Step 4)

Run in `code/04_metrics/`, in order:

1. `step4_1to3_subspace_metrics.R` — Grassmannian distance, principal angles, subspace/spectral recovery scores.
2. `step4_4to6_secondary_metrics.R` — trustworthiness, continuity, ARI, silhouette.
3. `step4_9_compile_results.R` — compiles the master results table.

**Estimated runtime:** ~48 minutes for step 1, several hours for step 2 (real-data phase runs serially due to memory constraints on large datasets).

---

## 8. Phase-space and failure-boundary analysis (Step 5)

Run in `code/05_phase_space/`, in order:

1. `step5_1_organize_by_axis.R`
2. `step5_2to7_axis_plots.R`
3. `step5_8_critical_thresholds.R`
4. `step5_9_phase_maps.R`

See `docs/step5_axis_findings.md` for the clipping-inertness and SymSim
batch-anomaly findings surfaced during this step's review.

---

## 9. Failure-mode detection (Step 6.1–6.8)

Run in `code/06_failure_modes/`, in this exact order (later scripts depend
on earlier ones' output tables):

1. `step6_1_technical_separation.R`
2. `step6_2_convert_null_control_replicates.R` → `step6_2_cluster_collapse.R`
3. `step6_3_discover_null_control_files.R` → `step6_3_phantom_clustering_v2.R`
4. `step6_4_variance_hijacking.R`
5. `step6_5a_pancreas_dpt.R`, `step6_5b_simulated_log2fc_discovery.R`, `step6_5c_simulated_variance_ratio.R`
6. `step6_6_prebuild_real_reference_spaces.R` → `step6_6_neighborhood_collapse.R`
7. `step6_7_subspace_rotation_slippage.R` → `step6_7b_slippage_gradient_analysis.R`
8. `step6_8_cross_reference_boundaries.R`

**Or, to run 6.3 through 6.8 unattended in one pass** (assumes 6.1/6.2 already
complete): see `code/06_failure_modes/run_step6_rerun.sh`.

**Estimated total runtime:** ~11–12 hours (Step 6.6 alone takes ~5–7 hours;
most other substeps are under 2 hours each).

## 10. Failure-mode retention review (Step 6.9)

No script — this is a synthesis/writing step. See
`docs/step6_9_failure_mode_review.md` for the final retained/caveated verdict
on all seven failure-mode categories.

---

## Debugging and provenance scripts

`code/dev_history/` contains every diagnostic, pilot-test, and rescue script
used while building and debugging the pipeline above. These are **not part
of the run sequence** — they are kept as provenance, documenting bugs found
and fixed along the way (see `PATCH_NOTES.md` for the narrative). Do not run
them as part of a fresh reproduction attempt.
