# Patch Notes

Technical changelog of bugs found and fixed across the project, one entry
per issue. For the full narrative context behind any entry, see
`PROJECT_HANDOVER.md`.

---

## SymSim — BiocGenerics S4 `rank()` conflict

**Affected file:** `vendor/SymSim_patched/R/simulation_functions.R`
**Upstream repo:** https://github.com/YosefLab/SymSim

### Problem
When BiocGenerics is loaded (via scran/scater), it registers an S4 generic
for `rank()`. SymSim's `Get_params()` calls `alply()` (from plyr) followed by
`do.call(c, temp)`, which produces an S4 list-like object. BiocGenerics' S4
method dispatch then validates this object before calling the rank method and
fails with:
    "error in evaluating the argument 'x' in selecting a method for function
    'rank': 'x' must be an array of at least two dimensions"
This error occurs on every `SimulateTrueCounts()` call regardless of
parameters because `Get_params()` is called unconditionally.

### Fix
Two changes to `R/simulation_functions.R`:
1. **`Get_params()` (~line 96):** Replace `alply(X, 1, ...)` + `do.call(c, temp)`
   with `as.numeric(t(as.matrix(X)))` — mathematically identical (row-major
   unrolling of matrix X) but produces a plain numeric vector that bypasses S4
   dispatch. Changed `rank(values)` to `base::rank(values)`.
2. **`SimulateTrueCounts()` (~line 666):** Changed
   `rank(rowSums(params[[3]][chosen_hge,]))` to
   `base::rank(as.numeric(rowSums(as.matrix(params[[3]][chosen_hge, , drop=FALSE]))))`.

### Reproduction
```r
Rscript R/01_simulation/install_symsim_patched.R
```
This script clones SymSim, applies both patches, and installs from source.
The `vendor/SymSim_patched/` directory is gitignored (231MB) but regenerated
by this script.

---

## Splatter — `mean.rate` does not control sparsity

**Affected files:** `R/01_simulation/param_dict.R`, `R/01_simulation/simulate_splatter.R`

### Problem
`mean.rate` was the initially assumed sparsity lever for Splatter. It has
essentially no effect — Splatter renormalizes gene means internally in a way
that cancels it out.

### Fix
`bcv.common` (biological coefficient of variation) identified as the real
lever via direct calibration sweep (`calibrate_splatter.R`), confirmed
monotonic against achieved sparsity. `param_dict.R`'s sparsity mapping keys
off `bcv.common`, not `mean.rate`.

---

## Splatter — dropout `low`/`high` labels inverted

**Affected files:** `R/01_simulation/param_dict.R`, `data/simulated/param_grid.csv`,
7,290 individual `.rds` files (Splatter main grid)

### Problem
`param_dict.R` mapped `"low" → dropout.mid=3.0` and `"high" → dropout.mid=1.0`.
A controlled sweep of `dropout.mid` (0.5 to 5.0, all else fixed) showed the
true relationship is the opposite: higher `dropout.mid` means more dropout
and less depth (achieved sparsity 0.9706 at mid=1.0 vs. 0.9925 at mid=3.0).
Every file labeled `dropout="low"` had genuinely more dropout than every
file labeled `"high"`. Symptom that led to the discovery: depth deviation
was 70.8% for `"low"` and only 40.0% for `"high"`, backwards from what the
labels implied.

### Fix
Since the underlying count matrices are valid — only the label was wrong —
files were relabeled, not regenerated. `param_grid.csv` backed up
(`param_grid_v1_INVALID_dropout_swap_backup.csv`), then `dropout` swapped
low↔high for the 3,645 rows in each direction. Each affected file's
`run_params$dropout` field patched to match (7,290 files, 0 errors).
`param_dict.R`'s `splatter_dropout` mapping corrected for future generation.

### Verification
Post-fix depth deviation ordering: `none` (0.74%) < `low` (40.0%) <
`high` (70.8%) — now monotonic in the expected direction.

---

## Splatter — depth calibration table was dropout-blind

**Affected files:** `R/01_simulation/calibrate_splatter_depth_dropout.R` (new),
`R/01_simulation/validate_output_inventory.R`

### Problem
`splatter_calib_depth.csv` only ever covered `dropout.type="none"`, despite
`calibrate_splatter.R`'s own design note acknowledging dropout layers on top
of sparsity. Depth validation was comparing `dropout=low/high` files against
a baseline that didn't apply to them.

### Fix
New script extends the original `lib.loc` sweep across all three dropout
levels using the corrected `dropout.mid` values, producing
`splatter_calib_depth_dropout.csv` (30 combinations: dropout_label × lib_loc
× actual_depth). `validate_output_inventory.R`'s `get_expected_depth()`
updated to use this table for Splatter.

### Verification
Flagged files dropped from 62% (6,811/10,935) to 37.3% (4,077/10,935); mean
deviation 37.2% → 17.3%. The `dropout=none` row of the new table exactly
reproduces the original table (consistency check passed).

---

## scDesign3 — `sparsity_label` and `dropout` had zero generative effect

**Affected files:** `R/01_simulation/simulate_scdesign3.R` (rewritten; original
preserved as `simulate_scdesign3_v1_INVALID_sparsity_inert.R`),
`R/01_simulation/simulate_null_control_scdesign3.R`

### Problem
The original script generated one count matrix per `fit_key` and saved it
unchanged across all 5 sparsity labels and all 3 dropout values sharing that
key. Confirmed via direct byte-for-byte comparison: 20/20 randomly sampled
fit_key groups showed identical `sum(counts)` across all 5 sparsity
replicates. Root cause: `dropout` mapped to `family_use` (nb/zinb switch),
but the calibrated `zero_inflation_pi` values were never actually passed
into the generation call; `sparsity_label` had no generation-time mechanism
at all.

### Fix attempted and rejected
`family_use="zinb"` refitting: `fit_marginal` (~1,027 sec) + `extract_para`
(~711 sec) ≈ 29 minutes for just two steps at the smallest test size
(projected ~5 days for the full 244-fit grid). The extracted zero-inflation
matrix returned entirely `NA` on this package version, crashing `simu_new()`
with a `dimnames` mismatch.

### Fix adopted
`family_use` fixed at `"nb"` for all fits. `dropout` and `sparsity_label`
reimplemented as calibrated post-hoc stochastic zero-masking (Bernoulli,
seeded per row), composed via `p_combined = 1-(1-p_dropout)(1-p_sparsity)`.
Sparsity ladder: `{0.7:0.00, 0.8:0.15, 0.9:0.35, 0.95:0.55, 0.98:0.75}`.
Dropout: `{none:0.00, low:0.10, high:0.40}`. `fit_key` simplified (dropout
no longer needs its own fit): unique fits dropped from 244 to 82.

### Verification
Zero degenerate empty cells tested at both grid extremes, including the
harshest combined setting (dropout=high + sparsity=0.98, p_combined=0.85).
Post-fix: rank-order broken 0/2,187 (was 2,187/2,187). Old data backed up
(6.0GB) before regeneration; full regeneration (82 fits, 10,940 files) ran
with 0 errors in 1h16m.

---

## scDesign3 / Splatter / SymSim null-control — missing `gene_strategy`/`clipping` metadata

**Affected files:** `R/01_simulation/simulate_null_control_scdesign3.R`,
`simulate_null_control_splatter.R`

### Problem
`data.frame()` construction crashed on `NULL` `gene_strategy`/`clipping`
fields while validating the 30 new null-control replicate files. Confirmed
as a generation-script gap, not a bad source value — `param_grid.csv` had
the correct values throughout. Splatter and scDesign3 were both affected;
SymSim's null-control script already had these fields correct.

### Fix
Patched at the script level (both affected scripts) so these fields are
populated from `param_grid.csv` as they should have been from the start.
10 files per affected simulator repaired.

---

## SymSim — `sparsity_label` had zero generative effect

**Affected file:** `R/01_simulation/simulate_symsim.R` (rewritten; original
preserved as `simulate_symsim_v1_INVALID_sparsity_inert.R`),
`R/01_simulation/simulate_null_control_symsim.R`

### Problem
The original script's own comments already documented this: "Sparsity
label does not affect simulation (Sigma fixed at 0.4)." Same underlying
mistake as the scDesign3 bug above — a parameter that seemed like it should
control sparsity didn't, and the axis was declared uncontrollable rather
than searched for the correct lever. `dropout` (via `alpha_mean`/
`depth_mean`) was already confirmed working and needed no change.

### Fix
Same post-hoc masking mechanism as scDesign3, applied only to
`sparsity_label`. `fit_key` unchanged (sparsity was never part of it, so
fit count stayed at 244).

### Verification
Zero empty cells tested at the worst-case setting (depth=500, dropout=high,
separability=high, smallest population size). Post-fix rank-order broken
0/2,187 (was 2,187/2,187). Depth deviation essentially unchanged (0.82%),
confirming the fix didn't disturb the already-correct dropout/depth
relationship. Old data backed up (22GB) before regeneration. The main-grid
regeneration was interrupted once (terminal application closed, killing the
process despite `nohup`) at 216/244 fits done; resumed cleanly from
checkpoint with zero data loss. The null-control regeneration for SymSim
was separately, accidentally missed during recovery from that same
interruption — caught by file-not-found warnings in validation, run
separately afterward.

---

## Validation script — `mclapply` fork instability with `Matrix`/`cholmod`

**Affected file:** `R/01_simulation/validate_output_inventory.R`

### Problem
Validation logic worked correctly in serial execution but crashed when run
under `mclapply` forking, specifically interacting badly with the `Matrix`
package's `cholmod` C library.

### Fix
Switched to plain `lapply` (serial). Full-grid validation (32,865 files)
completed in about 3 minutes serially, so parallelism wasn't providing a
meaningful speed benefit anyway.

---

## Splatter — 4 files traced to abandoned first-invocation RNG substream (run_ids 4, 5, 9, 3289)

**Affected files:** 4 `.rds` files in `data/simulated/splatter/`

### Problem
Ground-truth extraction (Step 1.8) found these 4 files' true-signal
subspace didn't match their paired production count matrices. Traced to
the very first invocation of `simulate_splatter.R`, run under
BiocParallel's `MulticoreParam` before it was replaced with
`parallel::mclapply` the same day due to a separate parallel-backend
defect. BiocParallel manages its own per-worker RNG substreams, which
plausibly overrode or interacted with the seed argument passed into
`splatSimulate()`. Confirmed isolated to this abandoned first invocation by
direct regeneration of 124 other files sharing the same early time window
(under the corrected engine) — all 124 matched correctly.

### Fix
Regenerated via the current `simulate_splatter.R` logic. Old files backed
up locally, not deleted.

### Verification
A third, independent `splatSimulate()` call reproduced the fix exactly for
all 4 run_ids.

---

## Splatter — 6 files with unresolved ground-truth mismatch (run_ids 3284, 3285, 6565, 6569, 9849, 9850)

**Affected files:** 6 `.rds` files in `data/simulated/splatter/` (excluded
from ground truth, not fixed)

### Problem
Found during the full 10,935-row ground-truth extraction run; confirmed
disjoint from the abandoned-invocation issue above. All 6 share sparsity in
{0.95, 0.98} and n_cells=200 (necessary but not individually sufficient —
1,458 rows share both conditions, only these 6 fail).

### Root cause
Not found. Every plausible explanation was checked directly and ruled out:
relabel-metadata inconsistency, dropout value/type, gross file corruption,
`gene_strategy` as a real causal input (proven invisible to
`splatSimulate()`), session/fork RNG carryover (ruled out both empirically
and mechanistically — Splatter v1.26.0 uses `withr::with_seed()`, confirmed
via source inspection), silent parameter drift in `param_dict.R` or
`simulate_splatter.R` (checked against git history, unchanged since first
commit), and package/R version drift (checked against `renv.lock` history,
completely static).

### Resolution
These 6 rows (0.055% of the main grid) were excluded from the ground-truth
set rather than shipped with a mismatched label — logged in
`splatter_unresolved.csv`, not silently dropped. Assessed as
disproportionate to chase further into Splatter's internal C/R
implementation given the confirmed, small, already-excluded scope (0.09%
of the full main grid checked at 100% coverage, no evidence of broader
extent).

---

## Step 3 preprocessing — 13,091-file corruption from disk-full + in-place overwrite

**Affected files:** ~13,091 `.rds` files across the 4 linear-family methods
(loadings backfill)

### Problem
`/mnt/extra` reached 100% capacity mid-write during the Step 3.1–3.4
loadings backfill, combined with an in-place-overwrite write pattern
(rather than atomic temp-file + rename), leaving torn files behind.

### Fix
All affected files regenerated. A later Step 3.10 full-grid integrity scan
found 2 additional files from this same incident (both real-data files —
Baron, Muraro) that had been missed in the original recovery scope; these
were separately fixed at that time.

### Verification
Full re-read of all 32,814 raw-PCA files confirmed 0 corruption after the
fix. Step 3.10's manifest scan (196,920 files) confirmed 0 unreadable files
project-wide.

---

## Step 3.10 — real-data Raw PCA `nv` mismatch and missing loadings

**Affected files:** 2 real-data `.rds` files (`pca_raw`)

### Problem
Found during the Step 3.10 full-grid integrity manifest build: two
real-data Raw PCA files used `nv=50` inconsistent with the rest of the
real-data grid, and were missing their loadings matrix.

### Fix
Regenerated with correct `nv` and loadings included. The fix is
self-documented in each file's own `flags`/`loadings_note` metadata field
(`"nv_mismatch_and_missing_loadings_fix_2026-08-20"`).

---

## SCTransform v2 — zero-total-count cells

**Affected files:** 16 `.rds` files (`pca_sctransform_v2`, simulated)

### Problem
A small number of cells in a small number of files had zero total UMI
count, which SCTransform's fitting cannot handle. Two rescue passes
(`scripts/step3_5_sctransform_rescue.R`, `_rescue2.R`) resolved most
affected files; 16 remained with a genuinely zero-count cell.

### Fix
The single affected cell explicitly excluded from that file's SCTransform
input. Exclusion logged with cell ID and reason directly in the file's own
output metadata, not silently dropped.

### Note
This produces a deliberate, expected one-cell size mismatch between these
16 files' embeddings and their source SCE (embedding has exactly 1 fewer
cell). This same pattern was independently rediscovered and re-confirmed
three times in Step 6 (as 17 rows, not 16 — see the Step 6.1 entry below)
via cross-file alignment checks in Steps 6.1, 6.4, and 6.6, each time
correctly identified as this same known, benign case rather than treated
as a new bug.

---

## PBMC68k — zero usable ground-truth labels

**Affected file:** `data/real/pbmc68k/pbmc68k_harmonized.rds`

### Problem
Going into Step 4, PBMC68k had 0/65,690 cells with usable ground-truth
cell-type labels (100% `NA` in `true_group`).

### Fix
Full re-annotation pipeline: Louvain clustering plus marker-based scoring,
with a CD3D-gating correction for a CD8⁺ T-cell/NK marker-overlap issue.
Production harmonized SCE file replaced; `data/real_data_inventory.csv`
metadata corrected; `embedding_manifest.csv` re-joined to propagate the fix.

### Verification
65,690/65,690 cells labeled across 5 canonical types (B_cell, CD4_T,
CD8_T, Monocyte, NK).

### Follow-on bug (see Step 6.1 entry below)
The `embedding_manifest.csv` re-join above was metadata-only — it never
touched the 6 individual pbmc68k embedding `.rds` files' own internally
stored `true_group` field, which remained 100% `NA` until caught and
patched separately in Step 6.1.

---

## scDesign3 — dropout label swap (Step 5 axis-sweep review)

**Affected file:** `scripts/step5_fix_scdesign3_dropout_labels.R`

### Problem
Visual review of Step 5's single-axis dropout recovery plot for scDesign3
showed an inconsistent pattern, traced to a dropout label swap distinct
from (and later than) the Step 1.7 scDesign3 sparsity-inert bug above.

### Fix
Corrected via the dedicated fix script; documented in
`docs/step5_axis_findings.md`.

---

## Step 6.1 — `aricode::ARI()`/`AMI()` silent character-to-integer miscoercion

**Affected files:** `scripts/step6_1_technical_separation.R` and every
subsequent Step 6 script using `aricode` (`step6_2`, `step6_3`, `step6_4`,
`step6_6`)

### Problem
`aricode::ARI()`/`AMI()` internally coerce label vectors via `as.integer()`
directly, not via `factor()`. Numeric-looking character strings (e.g.
`"1"`, `"2"`) survive this coercion; any other character labels (e.g.
`"CD4_T"`, `"Group1"`, `"alpha"`) become `NA`, producing:
    "Error in sort_pairs(c1, c2) : NA are not supported"
with an accompanying "NAs introduced by coercion" warning. Confirmed via a
direct test batch: only rows with purely numeric-string true_group values
(SymSim's `"1"`–`"5"` labels) succeeded; every other tested row failed.

### Fix
All label vectors passed through `as.integer(factor(x))` before any
`aricode::ARI()`/`AMI()` call, in every script that uses this package.

### Verification
Re-run of the same test batch: all 10/10 rows succeeded; recomputed ARI
matched Step 4.5's stored values in 6/10 cases exactly (using
`nstart=25`), confirming both the fix and the correct k-means
configuration to use going forward.

---

## Step 6.1 — `is_null_control %in% FALSE` silently drops all real-data rows

**Affected files:** `scripts/step6_1_technical_separation.R` and every
subsequent Step 6 script filtering on `is_null_control`

### Problem
`is_null_control` is `NA` (not applicable), not `FALSE`, for real-data
rows in both `embedding_manifest.csv` and `step4_master_results_table.csv`.
`%in%` does not propagate `NA` — `NA %in% FALSE` evaluates to `FALSE` — so
`manifest$is_null_control %in% FALSE` silently excluded 100% of real-data
rows (36/36) from every dry run without any error or warning.

### Fix
Changed to `!(manifest$is_null_control %in% TRUE)`, which correctly
includes both `FALSE` (simulated, non-null) and `NA` (real, not
applicable) rows while excluding only actual `TRUE` null-controls.

### Verification
Post-fix row count for `pbmc68k` and `tabula_sapiens_lung` in the
filtered working set: 6 rows each (as expected: 1 dataset × 6 methods),
vs. 0 before the fix.

---

## Step 6.1 — stale `true_group` in the 6 pbmc68k embedding files

**Affected files:** `data/processed/pca_{raw,libnorm,log,shiftedlog,
sctransform_v2,glmpca}/real/pbmc68k/*.rds` (6 files)

### Problem
All 6 pbmc68k embedding files had `true_group` = 100% `NA` (65,690/65,690),
even though the source SCE (`pbmc68k_harmonized.rds`) had correct,
fully-labeled data following the Step 4.7 re-annotation. Root cause: that
re-annotation's `embedding_manifest.csv` re-join was metadata-only (see the
PBMC68k entry above) and never touched these files' own internal
`true_group` field, which had been baked in at original file-creation time,
before the re-annotation existed.

### Fix
`true_group` patched in place in all 6 files (atomic write: temp file +
rename) with values pulled directly from the corrected source SCE,
positionally. The fix is logged into each file's own `flags` field. No
independent per-cell ID existed in either object to directly prove
positional alignment (neither has colnames/rownames on the relevant
matrices), so this was documented as an assumption, not a proven fact, at
the time of the fix.

### Verification
Recomputed k-means (`nstart=25`, seed 42) + ARI on the patched files
matched Step 4.5's stored ARI values in 5/6 methods to 4 decimal places
exactly, and the 6th (`pca_raw`) within 0.0004 — strong indirect
confirmation the positional-alignment assumption holds (a genuinely broken
alignment would produce near-random ARI regardless of the stored value).

---

## Step 6.3 — 30 null-control replicate files never preprocessed

**Affected files:** `scripts/step6_2_convert_null_control_replicates.R`
(new); no existing Step 3 script modified

### Problem
Only "replicate 1" of each null-control condition (the 5 rows already
embedded in `param_grid.csv` per simulator) had ever been run through
Step 3's 6 preprocessing methods. The 30 newer replicate files (Step 1.6,
`rep2`/`rep3`) lived in a separate raw-list format
(`data/simulated/null_control/`, distinct from `data/simulated/sce/`) that
Step 3's driver scripts never scanned.

### Fix
New standalone conversion script builds standardized SingleCellExperiment
objects from the 30 raw files, matching Step 1.8's schema, written into
`data/simulated/sce/{simulator}/` under new `_rep2`/`_rep3` filenames. All
6 existing, unmodified Step 3 scripts then simply re-run as-is — their
`skipped_existing` guards fast-skip the ~32,820 already-done files and
compute only the 30 new files × 6 methods = 180 new embeddings.

### Verification
30/30 output files confirmed present for all 6 methods, 0 errors.

---

## Step 6.3 — Calinski-Harabasz literal threshold (>10) is indiscriminate

**Affected file:** `scripts/step6_3_phantom_clustering.R` (v1, superseded
by `_v2.R`)

### Problem
The literal spec threshold flagged 270/270 (100%) null-control files —
suspicious given the phantom-clustering hypothesis was expected to be
method-dependent. A 20-draw pure-Gaussian-noise simulation (matched n=1000,
30 dims) showed CH scores of 26.6–27.8 even with zero real structure,
confirming CH has no built-in chance-correction and an absolute threshold
of 10 does not discriminate real structure from noise at this n/k.

### Fix
Empirical calibration: 200-draw noise baseline per embedding dimensionality
(pca_raw=50 dims, other 5 methods=30 dims), 95th percentile used as a
relative threshold instead of the literal absolute value.

### Verification
Under the corrected threshold, flag rate varies meaningfully by method
(64.4%–100%), unlike the literal threshold's uniform 100% — confirming the
calibrated version actually discriminates.

---

## Step 6.5a — `destiny`/`diffusionMap` installation failure

**Affected file:** `scripts/step6_5a_pancreas_dpt.R`

### Problem
`BiocManager::install("destiny")` failed with a cascading dependency chain:
`smoother` unavailable on CRAN for the installed R version, breaking `VIM`,
breaking `destiny` itself (`ERROR: dependencies 'VIM', 'smoother' are not
available for package 'destiny'`). `diffusionMap` was also checked and
found unavailable.

### Fix
Diffusion pseudotime implemented manually per Haghverdi et al. 2016
(Gaussian kernel affinity, symmetric-normalized Markov transition matrix,
`RSpectra::eigs_sym` for the leading eigenvectors, standard DPT distance
weighting `lambda/(1-lambda)`), using only already-available packages
(`RSpectra`, confirmed installed without issue).

---

## Step 6.5a — stale `pca_sctransform` directory accidentally globbed

**Affected file:** `scripts/step6_5a_pancreas_dpt.R`

### Problem
An initial recursive file-discovery glob picked up both
`data/processed/pca_sctransform_v2/` and a previously undocumented,
separate `data/processed/pca_sctransform/` (no `_v2`) directory, producing
21 rows instead of the expected 18. Directory timestamps confirmed
`pca_sctransform/` is ~3 weeks older and its files are half the size of
the current `_v2` output — a superseded, pre-loadings-backfill artifact
from an earlier Step 3.5 rescue attempt.

### Fix
File discovery restricted to the canonical 6 method directory names
explicitly, rather than a recursive glob.

### Note
The stale `data/processed/pca_sctransform/` directory itself was left in
place, not deleted, pending explicit confirmation — still an open item as
of this writing.

---

## Step 6.5c — circular ground-truth-comparison design (v1)

**Affected file:** `scripts/step6_5c_simulated_variance_ratio.R` (v1
discarded, replaced in place)

### Problem
v1 projected group-mean-transformed raw counts through each method's own
fitted loadings, then compared against the same method's actual embedding
centroids. This is mathematically circular: group-averaging and any fixed
linear projection commute unconditionally, so for any correctly-matched
linear transform the ratio is guaranteed to be ≈1.0 by matrix algebra, not
by measurement. Symptom: `pca_raw`, `pca_libnorm`, and `pca_log` all
returned ratios of 0.99999–1.00003 (6 decimal precision) — too perfect to
be a real measurement.

### Compounding bug in the same version
`pca_shiftedlog`'s transform was assumed to be `log1p(x+1)`. The actual
`step3_4_shiftedlog_pca_full.R` formula uses a data-driven, per-file shift
(`delta_opt` = median of nonzero library-normalized values), verified
directly from source. This wrong assumption produced `pca_shiftedlog`'s
anomalous 22.5% flag rate in v1 — a red herring that would have persisted
even if the circularity itself had gone unnoticed.

### Fix
Full redesign (not a parameter tweak): compares an ANOVA-style signal
fraction (between-group SS / total SS) computed independently in two
genuinely unrelated spaces — the method's transformed full gene-expression
matrix, and the method's actual embedding — with no shared projection
linking the two sides. The correct, verified `pca_shiftedlog` transform
formula (including the per-file `delta_opt`) used throughout.

### Verification
20-file pilot test: ratios spread meaningfully (0.55–14.7) instead of
clustering at exactly 1.0, confirming the circularity was removed.

---

## Step 6.6 — `apply()` on sparse matrix causes worker crashes

**Affected file:** `scripts/step6_6_neighborhood_collapse.R`,
`scripts/step6_6_prebuild_real_reference_spaces.R` (new, added as part of
this fix)

### Problem
HVG selection used `apply(logged, 1, var)` on a sparse `dgCMatrix`, which
has no sparse-aware path and forces expensive repeated dense conversions.
Caused 9/500 (later 4/500, after a partial fix) worker-process crashes in
`mclapply` — invisible to `tryCatch` since these were process deaths, not
R-level errors (`mclapply` warning: "N parallel function calls did not
delivered results"). Crashes concentrated on the two highest-gene-count
real datasets (Tabula Sapiens Lung: 56,139 genes; Tasic 2018: 42,865
genes) under concurrent 4-worker memory pressure.

### Fix
Two parts: (1) HVG variance computed via the sparse-efficient identity
`Var(X) = E[X²] - E[X]²` using `Matrix::rowMeans()` (properly sparse-aware),
replacing `apply()` entirely; (2) all 6 real datasets' reference spaces
pre-built sequentially, single-process, with explicit intermediate `gc()`
calls between steps, before the parallel run starts — removing the
concurrent-memory-pressure scenario for the highest-risk files entirely.

### Verification
Worker crashes: 9/500 → 4/500 (fix 1 alone) → 0/500 (fix 1 + fix 2
combined).

---

## Step 6.6 — cache-write race condition (non-atomic `saveRDS`)

**Affected file:** `scripts/step6_6_neighborhood_collapse.R`

### Problem
`build_reference_space()`'s cache write used plain `saveRDS()` with no
atomic temp-file-then-rename, violating this project's established
atomic-write convention (see the Step 3 corruption entry above). Multiple
workers processing different methods for the same underlying source file
could each observe `file.exists()==FALSE` simultaneously (a
time-of-check-to-time-of-use race) and write to the identical cache path
concurrently, producing a torn file. Symptom: 963 errors (0.49% of the
full run) with messages like "ReadItem: unknown type X, perhaps written
by later version of R" and "error reading from connection," clustering by
shared `run_id` across multiple methods simultaneously (e.g. run_id 3143:
`pca_raw`, `pca_sctransform_v2`, `pca_glmpca`, `pca_libnorm` all failed
together) — the signature of one corrupted shared file breaking every
subsequent reader, not independent random I/O noise.

### Fix
Cache write changed to temp-file + atomic rename (unique temp filename per
process, via PID). Cache read wrapped in `tryCatch`, falling through to a
full rebuild if the cached file is found to be corrupted on read.

### Resolution
A full integrity scan of the cache directory found 10/32,811 files
(0.03%) genuinely corrupted; these were deleted and their dependent rows
recomputed. Final error count after the fix and retry: 17, matching the
known-benign SCTransform zero-count-cell exclusion pattern (see that entry
above) exactly — confirmed via direct file-path set comparison against
Step 6.1's identical 17-row error set, not just matching counts.

---

## Step 6.7 — gene-count mismatch for SCTransform v2 / GLM-PCA

**Affected file:** `scripts/step6_7_subspace_rotation_slippage.R`

### Problem
167/500 (33.4%) dry-run errors, all "non-conformable arguments,"
concentrated exactly in `pca_glmpca` (82) and `pca_sctransform_v2` (85).
Initial hypothesis (the `gene_strategy` grid axis reducing gene count) was
tested directly and ruled out — an `hvg2000`-labeled `pca_raw` row had
full 10,000-gene loadings identical to `true_group_means`, confirming
`gene_strategy` is metadata-only for the 4 linear methods, no actual
subsetting occurs. Root cause: both `step3_5_sctransform_with_loadings.R`
and `step3_6_glmpca_simulated.R` have an internal `MAX_GENES=3000` cap
independent of `gene_strategy`, so for Splatter's 10,000-gene files these
two methods' `loadings` cover a smaller, method-specific HVG subset than
`true_group_means`'s full gene set.

### Fix
`true_group_means` subsetted to `loadings`' actual genes, matched by gene
name (not position), before computing the true subspace basis. Since
SCTransform v2 and GLM-PCA each select their own internal HVG subset
(potentially different from each other), the true-subspace cache is keyed
by `(source_file, method)`, not by `source_file` alone as it was for
Step 6.6's analogous cache.

### Verification
Dry-run error count: 167/500 → 0/500.

---

## Step 6.7b — uninformative "ceiling effect" diagnostic

**Affected file:** `scripts/step6_7b_slippage_gradient_analysis.R`

### Problem
An initial `always_above_30deg_ceiling_effect` flag fired for all 30/30
axis-method combinations tested, providing no discriminating information
between genuinely stressor-insensitive methods and genuinely
stressor-sensitive ones.

### Fix
Replaced with a range-of-movement metric (max angle − min angle across
each axis's stressor levels), which correctly separates the two cases
(e.g. `pca_raw`: range 0.97–5.3° across axes vs. `pca_shiftedlog`: range
up to 27.5° for separability).

---

## Repository — 32 stray root-level `.log` files

**Affected:** repository root (32 files removed), `.gitignore`

### Problem
32 raw verbose stdout redirects (e.g. `step3_6_run.log`,
`loadings_backfill_run.log`) had been committed directly at the repository
root rather than the tracked `logs/` directory, evidently from
`Rscript foo.R > run.log 2>&1`-style redirects run from the project root.
Total size was modest (204K) but they cluttered the repo root and had no
scientific value beyond what's already distilled into the tracked
`logs/step3_*_progress.csv` summaries.

### Fix
All 32 files removed via `git rm`. A root-level `*.log` ignore rule added
to `.gitignore` to prevent recurrence.

---

## Repository — stale `README.md`

**Affected file:** `README.md`

### Problem
Described the project structure as `R/01_simulation/` through
`R/04_analysis/` only — never mentioned `scripts/`, `docs/`, or any work
from Steps 3 through 6, despite those representing most of the project's
actual completed work.

### Fix
Rewritten to accurately describe the full current repository structure,
including both code-organization conventions (`R/` for Steps 1–2,
`scripts/` for Steps 3–6) and every tracked top-level directory.

---

## Repository — stale `PROJECT_HANDOVER.md`

**Affected file:** `PROJECT_HANDOVER.md`

### Problem
Document stopped at "EDA Checkpoint 1 — Complete (2026-07-24)," silently
missing Steps 2 through 6 entirely — roughly half the project's actual
completed work.

### Fix
Fully rewritten from scratch, covering Step 1.1 through Step 6.9 in full
narrative detail.

### Process note
The first replacement attempt (commit `b4d3f6a`) was a silent no-op: the
downloaded replacement file saved to the local `Downloads` folder under a
different filename than expected, so the `cp` command intended to
overwrite `PROJECT_HANDOVER.md` was never actually run, and git committed
the unchanged original content under a new commit hash. Caught by
checking the commit's insertion/deletion line counts (66 lines changed
across 3 files — far too small for a genuine full-document rewrite) rather
than assuming the commit succeeded. Corrected in the following commit
(`27ad579`) after verifying the actual file content and line count
matched expectation before committing.

---

## Step 6.2 conversion script — invalid `return()` inside non-function `tryCatch()`

**Affected file:** `scripts/step6_2_convert_null_control_replicates.R`

### Problem
The skip-logic used `return(list(status = "skipped_existing", ...))` inside
a `tryCatch({...})` expression block that is not itself inside a function
definition — invalid in R, throws:
    "no function to return from, jumping to top level"
Never triggered on the script's first run (no output files existed yet, so
the `file.exists()` branch was never taken). Triggered on every one of the
30 files the second time the script ran, once outputs already existed —
including 10 files that only appeared to "exist" because they were
NTFS-corrupted husks still listed by the filesystem after an unclean
shutdown (see the NTFS corruption entry below).

### Fix
Restructured with `if (...) {...} else {...}` so the `tryCatch`'s
expression block's last-evaluated value becomes the returned value
naturally, without an explicit `return()`.

### Verification
Re-run after the fix: 10 converted (the genuinely corrupted symsim files),
20 correctly skipped (the already-good splatter/scdesign3 files), 0 errors.

---

## `/mnt/extra` — NTFS metadata corruption, 25 files lost

**Affected:** 25 files under `data/processed/` (Step 6.3-6.8 output CSVs/RDS)
and `data/simulated/sce/symsim/` (10 null-control replicate SCE files)

### Problem
`du`/`cat` on specific files returned "Input/output error." SMART health
check on the underlying NVMe came back clean (overall-health PASSED, 0
media/data integrity errors, 4% lifetime wear used) -- ruling out failing
hardware. `/mnt/extra` was found mounted read-only; `journalctl` showed
the actual cause: `ntfs-3g[814]: Trying to read non-allocated mft records
(65781 > 65760): Illegal seek` -- NTFS Master File Table corruption. The
drive's own SMART log recorded 58 unsafe shutdowns over its lifetime,
consistent with this being caused by an unclean shutdown (NTFS-3G leaves
the filesystem read-only automatically as a protective measure once it
detects this class of corruption, rather than risk further damage).

### Fix
`sudo ntfsfix /dev/nvme0n1p4` after unmounting -- the standard Linux NTFS
metadata-repair tool. Output confirmed the underlying cause directly:
"Metadata kept in Windows cache, refused to mount" -- i.e. Windows Fast
Startup had left the volume in a hibernated state Linux's NTFS driver
correctly refused to touch until cleared. Repair completed successfully
(MFT/MFTMirr compared OK, journal cleared); filesystem remounted `rw`.

### Verification
A full, clean sweep of both `/mnt/extra` (65,700 files) and `/mnt/extra2`
(135,944 files) after the repair found zero further I/O errors anywhere,
confirming the corruption was confined to exactly the 25 files already
identified, not filesystem-wide. `ntfsfix` does not recover data that was
already unreadable before the repair -- all 25 files needed regeneration,
not just a remount. All 25 were downstream analysis outputs (Step 6.3-6.8
results, null-control replicate conversions), regenerable from
already-validated scripts; no raw simulated data, real data, or
ground-truth extraction files were affected.

### Regeneration
`run_step6_rerun.sh` (new, chains the 8 remaining Step 6 sub-analyses
sequentially, `nohup`+`disown` for terminal-independent execution) run
unattended, ~12 hours (03:20-15:22 IST). Zero new errors beyond the
already-documented 17 benign SCTransform exclusions (recurring identically
for the fourth time across Steps 6.1, 6.4, 6.6, and this rerun). Two
tasic2018 rows (pca_raw, pca_libnorm) were silently missing from the
initial Step 6.6 rerun output (a worker-process crash uncaught by
`tryCatch`, not a data problem -- confirmed via `embedding_manifest.csv`,
which showed both source files present and `status="ok"` upstream); fixed
via the script's existing `--retry_errors` mode on just those 2 rows
rather than a full 196,822-row recompute. All regenerated statistics
cross-checked line-by-line (not just aggregate rates) against the original
results across all 8 sub-analyses -- exact or near-exact matches
throughout (see the arithmetic-correction entry below for the one genuine
discrepancy found, which was in the original documentation, not the data).

**Recommended follow-up, not yet done:** disable Windows Fast Startup on
this machine (if dual-boot) to prevent recurrence. The ~20GB of orphaned
NTFS allocation this incident left behind (`df` still reports it as used;
`du` cannot walk it) requires a full Windows `chkdsk /f` to reclaim --
`ntfsfix` does not do this.

---

## `docs/step6_9_failure_mode_review.md` / `PROJECT_HANDOVER.md` — Step 6.5a arithmetic error

**Affected files:** `docs/step6_9_failure_mode_review.md`, `PROJECT_HANDOVER.md`

### Problem
Both documents stated Step 6.5a's real-data over-smoothing branch flagged
"9/18 (50%)" of dataset-method combinations. Discovered during the NTFS-
recovery regeneration above, when the regenerated `step6_5a_pancreas_dpt.csv`
was independently recounted directly (not narrated): muraro 6/6 TRUE +
baron 3/6 (raw, sctransform_v2, glmpca) + segerstolpe 1/6 (glmpca) = 10,
not 9. Every individual per-dataset, per-method claim in both documents
(GLM-PCA flagged in all 3 datasets; SCTransform v2 in 2/3; Baron/GLM-PCA
the single lowest correlation at rho~=0.177) was independently re-verified
against the data and confirmed correct -- isolating the error to a simple
counting mistake in the top-line summary figure, not a computation or
methodology error.

### Fix
Both documents corrected to 10/18 (55.6%), with a note in
`step6_9_failure_mode_review.md` documenting the correction and pointing
here.
