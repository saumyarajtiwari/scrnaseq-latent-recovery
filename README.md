# scRNA-seq Latent Subspace Recovery

Mapping latent biological subspace recovery in scRNA-seq preprocessing
under sparsity, depth, and dropout stress.

## Core question

Which scRNA-seq preprocessing methods preserve the true latent biological
subspace, and at what technical thresholds do they fail? This is a
latent-recovery study, not a benchmark: the central contribution is
mapping *when and why* preprocessing distorts biological structure, using
simulated (Splatter, scDesign3, SymSim) and real (PBMC 68k, Muraro, Baron,
Segerstolpe, Tabula Sapiens Lung, Tasic 2018) datasets across six
preprocessing/dimensionality-reduction methods.

## Repository structure

This repository spans two phases of the project with different code
organization conventions, kept as-is rather than retroactively unified,
since both contain validated, already-published-on work:

- **`R/`** — Steps 1–2 (simulation generation and calibration; real-data
  download, QC, and harmonization) and early exploratory analysis
  (`03_eda_checkpoint2/`). Organized into numbered subdirectories by
  project phase. Includes intentionally-preserved buggy historical
  versions (e.g. files suffixed `_v1_INVALID_sparsity_inert`) as
  documented evidence of bugs found and fixed, not dead code.
- **`scripts/`** — Steps 3–6 (preprocessing pipelines, subspace-recovery
  metrics, phase-space/failure-boundary analysis, and the seven-category
  failure-mode detection framework). Flat naming convention:
  `step{N}_{substep}_{description}.R`.
- **`docs/`** — Structured findings documents intended for direct reuse in
  the manuscript (e.g. `step5_axis_findings.md`,
  `step6_9_failure_mode_review.md`).
- **`data/`** — `real/`, `simulated/`, and `processed/` subdirectories.
  Raw and processed `.rds` data files are excluded from version control
  (see `.gitignore`) and live on local/external storage; a small number of
  critical, size-bounded artifacts (ground-truth extraction outputs, one
  documented bug-evidence sample) are tracked as explicit exceptions.
- **`results/tables/`** — Tracked CSV summary tables from exploratory data
  analysis (Step 1).
- **`figures/`** — Rendered phase-space heatmaps and other visual outputs
  (Step 5.9).
- **`logs/`** — Tracked per-file progress logs (`step3_*_progress.csv`,
  etc.) for the largest batch-processing runs, kept as a processing
  audit trail. Ad hoc verbose stdout redirects are excluded (`*.log`).
- **`renv.lock`** — Pinned R package versions for reproducibility (`renv`).
- **`PROJECT_HANDOVER.md`** — Structured, continuously-updated handover
  record of project state, decisions, and open items.
- **`PATCH_NOTES.md`** — Chronological log of fixes and corrections
  applied across the project.

## Methods

Six preprocessing/dimensionality-reduction methods are compared: Raw PCA,
Library-size-normalized PCA, Log-PCA, Shifted-Log-PCA,
Pearson-residual/SCTransform-style PCA, and GLM-PCA.

## Status

Steps 1 through 6 (through Step 6.9) are complete as of the latest commit.
See `PROJECT_HANDOVER.md` for full current status and next steps.
