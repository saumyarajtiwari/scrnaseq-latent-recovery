# scRNA-seq Latent Subspace Recovery

Mapping latent biological subspace recovery in scRNA-seq preprocessing
under sparsity, depth, and dropout stress.

## Core question

Which scRNA-seq preprocessing methods preserve the true latent biological
subspace, and at what technical thresholds do they fail? This is a
latent-recovery study, not a benchmark: the contribution is mapping *when
and why* preprocessing distorts biological structure, using simulated
(Splatter, scDesign3, SymSim) and real (PBMC 68k, Muraro, Baron,
Segerstolpe, Tabula Sapiens Lung, Tasic 2018) datasets across six
preprocessing/dimensionality-reduction methods, evaluated against seven
distinct failure-mode categories.

## How to run this

See [`RUNNING.md`](RUNNING.md) for the complete, ordered execution guide
from environment setup through the final failure-mode review.

## Repository structure

```
code/
├── 01_simulation/       Simulation generation and calibration (Steps 1.1-1.8)
├── 02_real_data/        Real-data download, QC, and harmonization (Step 2)
├── 03_preprocessing/    Six preprocessing/DR methods (Step 3)
├── 04_metrics/          Subspace-recovery metric computation (Step 4)
├── 05_phase_space/      Phase-space and failure-boundary analysis (Step 5)
├── 06_failure_modes/    Seven-category failure-mode detection (Step 6.1-6.8)
├── eda_checkpoint2/     Real-dataset technical characterization
├── eda_checkpoint3/     Representative embedding visualization
└── dev_history/         Preserved debug/pilot/rescue scripts (not part of
                         the run sequence -- kept as provenance; see
                         PATCH_NOTES.md for what each one was for)

data/
├── simulated/
│   ├── ground_truth/    True biological signal per simulated file (tracked)
│   ├── param_grid.csv           The invariant 10,940-row parameter grid
│   └── null_control_grid.csv    45-row null-control manifest
└── real_data_inventory.csv

results/
├── step1_eda/           EDA Checkpoint 1 outputs (see note below)
├── step2_eda/           EDA Checkpoint 2 outputs
├── step3_eda/           EDA Checkpoint 3 outputs (best/worst/null panels)
├── step4_metrics/       Subspace-recovery metric tables + embedding manifest
├── step5_phase_space/   Critical-threshold tables and phase-space figures
└── step6_failure_modes/ All nine failure-mode detection result tables

docs/
├── real_data_metadata_catalog.md
├── step5_axis_findings.md              Clipping-inertness / batch anomaly findings
├── step6_9_failure_mode_review.md      Final retention verdict, all 7 failure modes
└── eda_checkpoint1_reconstruction_note.md   Documents a known code gap (see below)

logs/                    Tracked per-run progress logs for the largest batch jobs
```

## A known, documented gap: EDA Checkpoint 1's script

EDA Checkpoint 1's outputs are tracked (`results/step1_eda/`), but the exact
code that produced them was run interactively and never saved as a
standalone script. This is explained fully, including how to reconstruct it
from already-verified logic reused elsewhere in the pipeline, in
[`docs/eda_checkpoint1_reconstruction_note.md`](docs/eda_checkpoint1_reconstruction_note.md).
No other pipeline stage has this gap.

## Methods

Six preprocessing/dimensionality-reduction methods are compared: Raw PCA,
Library-size-normalized PCA, Log-PCA, Shifted-Log-PCA,
Pearson-residual/SCTransform-v2-style PCA, and GLM-PCA. Seven failure-mode
categories are tested: Technical Separation, Cluster Collapse, Phantom
Clustering, Variance Hijacking, Over-Smoothing, Neighborhood Collapse, and
Subspace Rotation Slippage. See `docs/step6_9_failure_mode_review.md` for
the final, evidence-reviewed verdict on each.

## Reproducibility

- R package versions are pinned via `renv.lock`; run `renv::restore()` before anything else.
- Every stochastic script sets `set.seed(42)` at minimum; several use additional per-row deterministic seeds for reproducible-but-distinct random draws (documented in each script's header comments).
- **Known reproducibility gaps**, documented rather than hidden: SymSim's installed source is pinned to a GitHub branch HEAD, not a fixed commit SHA; EDA Checkpoint 2's `irlba()` calls were run without a preceding `set.seed()`.

## Data availability

Raw and processed `.rds` data files (~340GB total) are not tracked in this
repository — they live on local/external storage during active development.
This repository tracks the small, essential, and non-regenerable artifacts:
ground-truth extraction outputs, final result tables, and all code needed to
regenerate everything else from the pinned parameter grid and fixed seeds.

## License

Code in this repository is released under the MIT License (see `LICENSE`).

## Status

Steps 1 through 6 (through Step 6.9) are complete as of the latest commit.
See `PROJECT_HANDOVER.md` for full narrative detail and `PATCH_NOTES.md` for
a structured changelog of every bug found and fixed along the way.
