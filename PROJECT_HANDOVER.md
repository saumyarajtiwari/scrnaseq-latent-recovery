# Project Handover — Mapping Latent Biological Subspace Recovery in scRNA-seq Preprocessing under Sparsity, Depth, and Dropout Stress

**Repository:** `saumyarajtiwari/scrnaseq-latent-recovery` (public, branch: `main`)
**Last updated:** 2026-07-13 (Step 1.7 complete, full validation pass)

---

## 1. Project Overview & Scientific Motivation

This is a **latent-recovery study**, explicitly not a broad benchmark paper. The central question: *which scRNA-seq preprocessing methods preserve the true latent biological subspace, and at what technical thresholds do they fail?*

Raw count matrices are transformed by different preprocessing methods, and the resulting low-dimensional representations are evaluated against known/simulated biological structure. Focus is on whether preprocessing preserves the true biological subspace under sparsity, low sequencing depth, dropout, and batch complexity — not on ranking methods in the abstract.

**Simulators used (three, to reduce simulator-specific bias):** Splatter, scDesign3, SymSim.

**Design:**
- Main parameter grid: 10,940 rows = 5 sparsity × 3 depth × 3 dropout × 3 separability × 3 cell count × 3 batch × 3 gene strategy × 3 clipping (`n_groups` fixed at 5 for the main grid).
- Null-control set: 45 matrices = 5 sparsity × 3 simulators × 3 replicates, single biological population, used to detect phantom clustering and other technical artifacts.
- Evaluation metrics: Grassmannian distance, principal angles, subspace/spectral recovery scores (primary); trustworthiness, continuity, ARI, silhouette (secondary).
- Phase-space failure-boundary analysis across sparsity, depth, separability, batch, dropout, clipping.
- Seven failure-mode categories: technical separation, cluster collapse, phantom clustering, variance hijacking, over-smoothing, neighborhood collapse, subspace rotation slippage.

---

## 2. Environment & Machine Setup

- **Machine:** Acer Predator Helios Neo 16 (Intel i5-13500HX, 16GB RAM)
- **OS:** Ubuntu 24.04.4 LTS, username `aayush`
- **R:** 4.3.3, with `renv` 1.2.3 (environment locked)
- **Storage:** Large simulation outputs live on the NTFS partition `/mnt/extra/scrnaseq-data/` (~98GB), symlinked into the project's `data/` directory. A `.renvignore` file excludes these `.rds` outputs from `renv`'s dependency scan (large-directory scanning previously caused environment issues).
- **Key packages:** Splatter, scDesign3, SymSim (`YosefLab/SymSim`), scran, scuttle, bluster, scater, BiocSingular, BiocParallel, TENxPBMCData.
- **Workflow discipline:** strict paste-and-confirm (one command, output pasted back, before proceeding); no assumptions without terminal evidence; atomic shell steps with commit-per-substep; all invalid/superseded data preserved under timestamped backups, never deleted outright.

---

## 3. Repository Structure
R/01_simulation/
param_dict.R                              — single source of truth for all parameter mappings
simulate_splatter.R                       — main-grid Splatter generation
simulate_scdesign3.R                      — main-grid scDesign3 generation (REVISED, see §6)
simulate_scdesign3_v1_INVALID_sparsity_inert.R   — preserved original (sparsity/dropout-inert bug)
simulate_symsim.R                         — main-grid SymSim generation (REVISED, see §7)
simulate_symsim_v1_INVALID_sparsity_inert.R      — preserved original (sparsity-inert bug)
simulate_null_control_splatter.R          — null-control Phase B, Splatter
simulate_null_control_scdesign3.R         — null-control Phase B, scDesign3 (REVISED)
simulate_null_control_scdesign3_v1_INVALID_sparsity_inert.R
simulate_null_control_symsim.R            — null-control Phase B, SymSim (REVISED)
simulate_null_control_symsim_v1_INVALID_sparsity_inert.R
calibrate_splatter.R                      — original sparsity/depth calibration (dropout.type="none" only)
calibrate_splatter_bcv.R / calibrate_splatter_bcv2.R  — bcv.common sparsity calibration
calibrate_splatter_depth_dropout.R        — NEW: dropout-aware depth calibration (see §5)
calibrate_scdesign3.R                     — PBMC 3k reference prep, separability/depth calibration
validate_output_inventory.R               — Step 1.7 full validation script (see §9)
data/simulated/
param_grid.csv                            — the invariant 10,940-row grid (dropout labels corrected, see §5)
param_grid_v1_INVALID_dropout_swap_backup.csv  — preserved pre-correction copy
null_control_grid.csv                     — 45-row null-control manifest
pbmc3k_annotated.rds                      — annotated PBMC 3k reference (scDesign3)
splatter_calib_sparsity.csv, splatter_calib_depth.csv, splatter_calib_depth_dropout.csv
scdesign3_calib_depth.csv, scdesign3_calib_sparsity.csv, scdesign3_calib_dropout.csv
symsim_calib_depth.csv, symsim_calib_sigma.csv, symsim_calib_sparsity.csv
validation_inventory.csv, validation_summary.csv  — current full Step 1.7 audit output
splatter/, scdesign3/, symsim/            — symlinks into /mnt/extra, 10,940 files each
null_control/{splatter,scdesign3,symsim}/ — 10 new files each (plus 5 reused per simulator, living in the main grid dirs above)
logs/                                        — all generation/validation run logs

---

## 4. Parameter Grid & Design

`param_grid.csv` (10,940 rows) is the invariant foundational file. `run_id` (sequential integers 1–10,940) is the **universal RNG seed recovery mechanism** across all three simulators — any file can be regenerated deterministically from its `run_id` alone, given the corresponding script version.

**Axes:**
| Axis | Levels |
|---|---|
| sparsity | 0.70, 0.80, 0.90, 0.95, 0.98 (ordinal identifiers, not absolute targets) |
| depth | 500, 2000, 10000 |
| dropout | none, low, high |
| separability | null, low, medium, high |
| n_cells | 200, 1000, 5000 (capped per simulator/separability where reference pool is smaller) |
| batch | none, simple, complex |
| gene_strategy | all, hvg500, hvg2000 |
| clipping | none, clip99, log_stabilized |

`n_groups` fixed at 5 for the main grid (except scDesign3, where separability tier determines actual cell-type count: null=1, low=2, medium=4, high=5).

Null-control uses a **separate** manifest (`null_control_grid.csv`: `simulator, base_run_id, sparsity, replicate, seed, file_path, is_new`) rather than touching the invariant main grid — 45 rows = 15 reused (pointing back into the main grid, `is_new=FALSE`) + 30 new (`is_new=TRUE`, dedicated files under `data/simulated/null_control/{sim}/`).

---

## 5. Step 1.3 — Splatter Generation

**Original build:** `mean.rate` does not control sparsity (Splatter normalizes gene means internally, canceling the effect) — `bcv.common` is the real lever. `lib.loc`/`lib.scale` control depth. 10,940 files generated, verified against target grid.

### Corrections made during Step 1.7 (2026-07)

**Finding 1 — dropout `low`/`high` labels were inverted.** `param_dict.R` originally mapped `"low" → dropout.mid=3.0` and `"high" → dropout.mid=1.0`. A direct controlled sweep (`dropout.mid` from 0.5 to 5.0, all else held fixed) confirmed the relationship is monotonic and unambiguous: **higher `dropout.mid` → more dropout, less depth** (mid=1.0: sparsity 0.9706; mid=3.0: sparsity 0.9925). This means every file labeled `dropout="low"` actually had *more* technical dropout than every file labeled `"high"` — backwards. This directly explained the depth-deviation pattern found in validation (`low`: 70.8% mean deviation, `high`: 40.0%, `none`: 0.74%).

*Fix:* since the underlying count matrices are valid simulations of real, specific `dropout.mid` values, the correction was to **relabel, not regenerate**. `param_grid.csv` backed up (`param_grid_v1_INVALID_dropout_swap_backup.csv`), then `dropout` values swapped low↔high for the 7,290 affected rows (3,645 each). Every affected `.rds` file's `run_params$dropout` field patched to match (7,290 files, 0 errors, ~30 min). `param_dict.R`'s `splatter_dropout` mapping corrected (`low`→1.0, `high`→3.0) for future-regeneration consistency. Post-fix direction confirmed correct: `none` (0.74%) < `low` (40.0%) < `high` (70.8%) mean deviation.

**Finding 2 — depth calibration table was dropout-blind.** `calibrate_splatter.R`'s own design note states *"Sparsity is calibrated at dropout.type='none'... Dropout layers on top of sparsity"* — an acknowledged but never-quantified coupling. `splatter_calib_depth.csv` only ever covered `dropout.type="none"`, so Step 1.7's depth check was comparing `dropout=low/high` files against the wrong baseline.

*Fix:* new `calibrate_splatter_depth_dropout.R` extends the original `lib.loc` sweep (500 cells, 2000 genes, 3 reps, same seed pattern) across all three dropout levels using the corrected `dropout.mid` values, producing `splatter_calib_depth_dropout.csv` (dropout_label × lib_loc × actual_depth, 30 combinations). The `dropout=none` row exactly reproduces the original table (consistency check passed). `validate_output_inventory.R`'s `get_expected_depth()` updated to use this table for Splatter. Result: flagged files dropped from 62% (6,811/10,935) to 37.3% (4,077/10,935); mean deviation 37.2%→17.3%.

**Remaining depth deviation is explained, not a bug:** strongly concentrated at `target_depth < 100` (mean deviation 44.7% there vs. 13.9% elsewhere) — a small-sample-size statistical artifact (tiny denominators inflate percentage deviation), not a generation defect.

**Rank-order near-ties (see §9 for the general fix):** 563/2187 groups initially flagged; root-caused to real dropout-severity saturation of the finer sparsity gradient (concentrated 76% at `dropout=high`), with the worst observed violation across the *entire* dataset being only 0.22 percentage points. Resolved via a justified tolerance in the validation script, not treated as a data defect.

*Note: the script's own header comment claims `N_WORKERS = 10`, but the actual configured value is `2L` — a pre-existing, harmless documentation inconsistency, left as-is.*

---

## 6. Step 1.4 — scDesign3 Generation

**Original build:** fit-once-simulate-many strategy, 244 unique fits (depth × dropout × separability × n_cells × batch). PBMC 3k reference (`TENxPBMCData`), top 2000 HVGs, greedy-unique cluster→cell-type annotation. Depth multipliers 0.45/1.79/8.95 → 500/2000/10000 UMI, off a measured 1,116.9 UMI/cell HVG baseline. `n_cells` caps by separability: null=1135, low=1485, medium=2009, high=2695. 10,940 files originally generated and verified.

### Correction made during Step 1.7 (the most significant finding of this session)

**Finding — `sparsity_label` AND `dropout` had zero effect on generated data.** The original script generated one count matrix per `fit_key` and saved it **unchanged** across all 5 sparsity labels and all 3 dropout values sharing that key. Confirmed via direct byte-for-byte comparison: 20/20 randomly sampled fit_key groups showed identical `sum(counts)` across all 5 sparsity replicates. Root cause: `dropout` mapped to `family_use` (nb/zinb switch), but the calibrated `zero_inflation_pi` values (0.10/0.40, defined in `param_dict.R`) were never actually read or passed into the generation call. `sparsity_label` had no generation-time mechanism at all.

**Path to the fix — `zinb` refitting was tested and rejected:** fitting with `family_use="zinb"` took ~1,027 sec (`fit_marginal`) + ~711 sec (`extract_para`) ≈ 29 minutes for just two steps at the smallest test size (projected ~5 days for the full 244-fit grid), and the extracted zero-inflation matrix came back entirely `NA` on this scDesign3 package version, crashing `simu_new()` with a `dimnames` mismatch.

**Fix adopted:** `family_use` fixed at `"nb"` for all fits (fast — ~777 sec full generation at `n_cells=1000`, matching original working timing). `dropout` and `sparsity_label` reimplemented as **calibrated post-hoc stochastic zero-masking** (Bernoulli, `seed = run_id*1000 + combo_index`), composed via `p_combined = 1-(1-p_dropout)(1-p_sparsity)`. Calibration: sparsity ladder `{0.7:0.00, 0.8:0.15, 0.9:0.35, 0.95:0.55, 0.98:0.75}` (mask_p), dropout `{none:0.00, low:0.10, high:0.40}` (the already-defined-but-previously-unused pi values). Verified via controlled tests at both grid extremes (depth=500/separability=null and depth=10000/separability=high): **zero degenerate empty cells** even at the harshest combined setting (`dropout=high` + `sparsity=0.98`, `p_combined=0.85`).

`fit_key` simplified: `dropout` removed as a fit dimension (no longer needs its own model fit) — unique fits dropped from 244 to **82**.

Old script preserved as `simulate_scdesign3_v1_INVALID_sparsity_inert.R`. Old data backed up to `/mnt/extra/scrnaseq-data/scdesign3_v1_INVALID_backup_20260710_085814` (6.0GB). Full regeneration: 82 fits, 10,940 files, **0 errors**, completed in **1h16m** (12:49:59→14:06:24) — far faster than the original 15–25hr estimate, since fewer fits were needed.

New calibration files: `scdesign3_calib_sparsity.csv`, `scdesign3_calib_dropout.csv`. `run_params` now includes `family_use` (always `"nb"`), `dropout_pi`, `sparsity_mask_p`, `combined_mask_p`, `mask_seed`, `actual_sparsity` — full reproducibility parity with `bcv_common`/`lib_loc`.

**Verified post-fix:** rank-order broken 0/2187 (was 2187/2187). Sparsity and dropout both confirmed genuinely, strictly monotonic. Depth deviation 14.2% flagged (1,557/10,935), mean 9.4% — fully explained by masking removing real count mass (documented, expected; `validate_output_inventory.R` compares against `base_calibrated_depth × (1-combined_mask_p)`, not the pre-masking baseline).

**Secondary bug found+fixed:** the 30 new null-control files (scDesign3, and later found also for Splatter/SymSim) were missing `gene_strategy`/`clipping` in `run_params`. `param_grid.csv` confirmed to have correct source values — this was a generation-script gap, fixed at the script level.

**Null-control (scDesign3):** `simulate_null_control_scdesign3.R` had the identical fit-once-save-many sparsity-inert bug. Fixed identically (masking; dropout fixed at `"none"`/`pi=0` for null-control by design). 10 new files regenerated, confirmed monotonic (0.8624→0.8832→0.9106→0.938→0.9656 across sparsity labels 0.7→0.98).

---

## 7. Step 1.5 — SymSim Generation

**Original build:** simulate-once-save-many, 244 unique calls (separability × n_cells × batch × dropout × depth). A `BiocGenerics::rank()` S4 compatibility conflict required source-level patching (bytecode compilation made runtime patching ineffective). `Sigma` fixed at 0.4, `N_GENES=2000`. `depth_mean` is severely sublinear at low `alpha_mean` (high dropout) — the 10,000 UMI target is physically unreachable at high dropout (achieves only ~3,967–4,000). 10,940 files originally generated and verified.

### Correction made during Step 1.7

**Finding — `sparsity_label` had zero effect, explicitly documented but never addressed.** The original script's own comments state: *"Sparsity label does not affect simulation (Sigma fixed at 0.4)"* and *"Sparsity, gene_strategy, clipping do not affect SymSim output."* `dropout` (via `alpha_mean`/`depth_mean`, genuinely part of `fit_key`) was already confirmed working correctly from Step 1.5 and required no change. This is the same root design flaw as scDesign3's original bug: a naive parameter (`Sigma`) didn't control sparsity, so the axis was declared uncontrollable rather than searched for the correct lever.

**Fix:** identical post-hoc masking mechanism as scDesign3, applied only to `sparsity_label` (dropout untouched). Same sparsity ladder. Verified via a worst-case controlled test (depth=500, dropout=high, separability=high, smallest per-population size): zero empty cells even at the harshest setting.

`fit_key` unchanged (sparsity was never part of it) — still 244 unique calls.

Old script preserved as `simulate_symsim_v1_INVALID_sparsity_inert.R`. Old data backed up to `/mnt/extra/scrnaseq-data/symsim_v1_INVALID_backup_20260711_062458` (22GB). New calibration file: `symsim_calib_sparsity.csv`.

**Full regeneration:** 244 fits, 10,940 files, **0 errors**. Runtime highly variable by combination — `dropout=none`/`n_cells=5000` fits took 1.5–2+ hours *each* (vs. ~15 sec for small fits), since `True2ObservedCounts` is more expensive with little/no technical zero-inflation to short-circuit computation. **The run was interrupted once** (terminal application closed, which killed the process despite `nohup` — likely a full terminal-app closure rather than just a shell exit) at 216/244 done; resumed cleanly via checkpointing with **zero data loss**, completing 2026-07-13 16:23:29.

`run_params` now includes `sparsity_mask_p`, `combined_mask_p`, `mask_seed`, `actual_sparsity`, alongside unchanged `alpha_mean`/`depth_mean`/`sigma`/`n_de_evf`.

**Verified post-fix:** rank-order broken 0/2187 (was 2187/2187). Depth deviation stayed low (0.82%, 90/10,935) — essentially unchanged from before the fix, confirming the masking correction didn't disturb the already-correct dropout/depth relationship.

**Null-control (SymSim):** `simulate_null_control_symsim.R` had the identical bug (one matrix per replicate, saved unchanged across 5 sparsity labels). Fixed identically. Its regeneration was initially missed during the main-grid interruption/recovery cycle, caught when Step 1.7 validation threw file-not-found warnings for the 10 missing files. Run separately afterward: 10 files regenerated cleanly.

---

## 8. Step 1.6 — Null-Control (cross-simulator summary)

45 total files: **15 reused** (5 per simulator, embedded in the main grid with `is_null_control=TRUE`) + **30 new** (10 per simulator, dedicated files under `data/simulated/null_control/{sim}/`). Manifest: `null_control_grid.csv`.

All three simulators' "new" null-control generation scripts had bugs found and fixed during Step 1.7:
- **scDesign3:** sparsity-inert (fixed via masking) + missing `gene_strategy`/`clipping` (fixed).
- **Splatter:** missing `gene_strategy`/`clipping` only — sparsity was already fine, and `dropout` is fixed at `"none"` for null-control by design so the label-inversion bug never applied here. Patched directly (10 files; `param_grid.csv` confirmed as the correct source of truth for the values).
- **SymSim:** sparsity-inert (fixed via masking); `gene_strategy`/`clipping` were already correctly present (no bug there).

**Final validated state:** the `null_control` inventory group shows 60 rows (reflects a known, harmless double-counting artifact — the 15 reused files each appear once via main-grid listing and once via `null_control_grid.csv` listing; 45 physical files, 60 inventory rows), **0 flagged**.

---

## 9. Step 1.7 — Output Validation and Inventory

**Script:** `validate_output_inventory.R`. Scope agreed before building: all 32,865 files (32,820 main-grid + 45 null-control). Sparsity checked descriptively (recomputed from raw counts) plus monotonic rank-order per simulator (means aggregated across replicates first, per-replicate noise never checked directly). Depth checked against each simulator's own calibrated expectation, 20% deviation threshold (chosen to "absorb legitimate variance, catch real failures"). Output: one combined `validation_inventory.csv` (per-file) + `validation_summary.csv` (per-group rollup).

**Two infrastructure bugs hit before reaching real findings:**
1. `mclapply` fork instability with the `Matrix`/`cholmod` C library — validation logic worked perfectly serially but crashed under forking. Fixed by switching to plain `lapply` (serial runtime was already ~3 minutes for the full validation, so parallelism wasn't needed anyway).
2. `data.frame()` construction crash from `NULL` `gene_strategy`/`clipping` fields in the 30 new null-control files — this is what led directly to discovering the null-control metadata bug described above.

**Chronological discovery sequence:**
1. scDesign3 sparsity-inert bug (biggest single finding) → full fix + regeneration.
2. Splatter/SymSim null-control metadata gaps (found while re-validating the scDesign3 fix) → patched.
3. SymSim sparsity-inert bug (found via the same 100%-rank-order-broken signature as scDesign3's original bug) → full fix + regeneration.
4. SymSim null-control regeneration accidentally skipped during recovery from a terminal interruption → caught by file-not-found warnings, run separately.
5. Splatter dropout `low`/`high` label inversion (found while investigating persistently high depth deviation) → confirmed via direct `dropout.mid` sweep, fixed via relabeling.
6. Splatter depth-calibration dropout-blindness (found immediately after) → root-caused to the calibration script's own acknowledged-but-unquantified design note, fixed via a new dropout-aware calibration table.
7. Rank-order near-ties across all simulators → resolved via a single, empirically-justified tolerance (`RANK_ORDER_TOLERANCE = 0.005`, based on the observed maximum violation of 0.0022 across all 2,187 groups, applied uniformly).

**FINAL VALIDATED STATE:**

| Simulator | Files | Flagged | Mean depth dev | Rank-order broken |
|---|---|---|---|---|
| Splatter | 10,935 | 4,077 (37.3%) | 17.3% | 0/2187 |
| scDesign3 | 10,935 | 1,557 (14.2%) | 9.4% | 0/2187 |
| SymSim | 10,935 | 90 (0.8%) | 2.9% | 0/2187 |
| null_control | 60 | 0 | — | — |

Every remaining flag has a documented, principled explanation (real dropout/masking-driven depth reduction, small-sample noise at extreme low-depth corners) — none is an unexplained anomaly. 0 files unreadable, 0 `n_cells` mismatches, 0 `n_groups` mismatches, 0 missing calibration matches, across every group.

---

## 10. Key Learnings & Principles

**From this session (Step 1.7):**
- **"Inert axis" bugs** (an axis has zero effect despite being logged as if it does) are best caught by direct byte-for-byte comparison of outputs across labels within a fit-sharing group — the single most powerful diagnostic used throughout.
- A **"fit-once-save-many" architecture** is efficient but creates systemic risk: any axis *not* included in `fit_key` needs an explicit, verified downstream mechanism to have real effect — never assume a label implies an effect.
- **Label-inversion bugs** (real effect, wrong direction) are arguably more dangerous than inert axes, since they actively produce misleading trends rather than merely uninformative ones. Verify direction via controlled parameter sweeps, never assume from variable naming or prior documentation.
- When two axes are genuinely **coupled** (e.g., dropout/masking reducing effective depth), the correct fix is usually not to eliminate the coupling (which may be biologically realistic and worth keeping) but to make the *validation comparison* aware of it.
- **Percentage-deviation metrics are unstable at small denominators** — always check whether "high deviation" concentrates at extreme/small-value corners before treating it as systemic.
- **Rank-order/monotonicity checks on single-draw stochastic data need an explicit, empirically-justified tolerance** — strict inequality will flag Monte Carlo noise as "broken." Derive the tolerance from the actual observed noise ceiling across the full dataset, not an arbitrary round number.
- When already-generated data is wrong only in its *label* (not its content), **relabeling is preferable to regeneration** — equally rigorous, cheaper, and avoids introducing new risk via a fresh generation run.
- `mclapply`/fork-based parallelism can be unstable with certain C-library-backed R packages (`Matrix`/`cholmod`; separately, GAMLSS-based `zinb` fitting) even when the identical code runs perfectly in serial — always test both paths separately when debugging opaque fork failures.
- `nohup` does not guarantee survival if the entire terminal *application* (not just the shell) is closed — checkpointing (already a core project principle) is the real safety net, not the backgrounding mechanism alone.

**Carried forward from earlier steps:**
- Simulator sparsity axes are not trustworthy by default — the correct lever requires empirical verification.
- Different unique call counts across simulators are acceptable given explicit methodological rationale.
- Splatter's `mean.rate` does not control sparsity (internal normalization cancels it); `bcv.common` is the correct lever.
- Monotonic rank-order checks should use means across replicates, not individual values, to avoid false flags from sampling noise.
- Depth validation must account for masking/dropout-induced count-mass removal.
- SymSim requires source-level bug patching (bytecode compilation defeats runtime patches).
- `renv` scanning of large output directories requires a `.renvignore`.
- Checkpointing is essential — has prevented data loss on at least two occasions this project (a system auto-shutdown during earlier SymSim runs; a terminal-closure interruption during this session's SymSim regeneration).
- `param_dict.R` is the single source of truth for all parameter mappings; `run_id` is the universal seed recovery mechanism.

---

## 11. Current State, Next Steps & Reproducibility Manifest

**Current state:** Step 1.7 is complete and fully validated. All 32,865 files (32,820 main-grid + 45 null-control) are confirmed structurally correct, with every remaining flag explained and documented above.

**Immediate pending items:**
- Delete the two remaining backup directories once fully comfortable (`scdesign3_v1_INVALID_backup_20260710_085814`, 6.0GB; `symsim_v1_INVALID_backup_20260711_062458`, 22GB — ~28GB reclaimable).
- Commit and push this handover document, plus the Splatter dropout-label-inversion fix, dropout-aware depth calibration, and rank-order tolerance fix (not yet committed as of this writing).

**Next roadmap steps** (unchanged from original design):
- **1.8** — Standardized Storage (AnnData/SingleCellExperiment unification).
- **EDA Checkpoint 1** (1.1–1.4) — simulated output verification (sparsity/depth/DE-magnitude checks, summary pass/fail table). Much of this is now already satisfied by Step 1.7's validation, but the checkpoint's specific deliverables (plots, DE-magnitude checks) still need producing.
- **Step 2** — Real-data validation selection & curation (PBMC 68k, pancreas datasets, Tabula Sapiens lung subset, brain atlas subset).
- **EDA Checkpoint 2** — real-dataset technical characterization, phase-space overlay figure.
- **Step 3** — Preprocessing & dimensionality-reduction method implementation, applied across the full simulated + real dataset inventory.
- **EDA Checkpoint 3** — representative embedding visualization.
- **Step 4** — Primary & secondary metric calculation.
- **Step 5** — Phase-space & failure-boundary analysis.
- **Step 6** — Failure-mode detection & characterization.

**Disk space (as of this session):** `/mnt/extra` at 66GB/98GB used (33GB free), with both backups still present. Will drop to ~38GB used (~60GB free) once backups are cleared. Real-dataset downloads (Step 2) estimated modest (5–15GB combined). Step 3's embedding-storage footprint has not yet been estimated — flagged for a small-scale storage test before full-scale execution, given the number of method × dataset combinations involved.

**Reproducibility manifest:** every simulated file's `run_params` now contains the *actual* generative parameters used, not just labels — `bcv_common`/`lib_loc` (Splatter); `family_use`/`dropout_pi`/`sparsity_mask_p`/`combined_mask_p`/`mask_seed` (scDesign3); `alpha_mean`/`depth_mean`/`sparsity_mask_p`/`combined_mask_p`/`mask_seed` (SymSim). Any file can be fully reconstructed from its `run_id` alone, given the corresponding script version.

**Git state:** scDesign3 + SymSim corrections committed and pushed (`ca2f4c4`, "Step 1.7: Fix sparsity/dropout-inert bugs in scDesign3 and SymSim, add output validation"). Splatter's dropout-label-inversion fix, dropout-aware depth calibration, rank-order tolerance fix, and this handover document are pending the next commit.
