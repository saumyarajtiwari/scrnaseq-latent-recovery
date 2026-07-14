# =============================================================================
# param_dict.R
# Simulator-specific numerical mappings for all qualitative grid labels.
# Sourced by all simulation scripts. Not a standalone executable.
# =============================================================================

# -----------------------------------------------------------------------------
# SPLATTER
# -----------------------------------------------------------------------------

splatter_separability <- list(
  "null"   = list(nGroups = 1, de.prob = 0.00, de.facLoc = 0.0, de.facScale = 0.0),
  "low"    = list(nGroups = 5, de.prob = 0.05, de.facLoc = 0.2, de.facScale = 0.4),
  "medium" = list(nGroups = 5, de.prob = 0.20, de.facLoc = 0.5, de.facScale = 0.4),
  "high"   = list(nGroups = 5, de.prob = 0.50, de.facLoc = 1.0, de.facScale = 0.4)
)

# CORRECTED 2026-07 (Step 1.7 audit): low/high were inverted. Direct
# empirical test confirmed dropout.mid is monotonic with dropout severity
# (higher mid -> more sparsity, less depth): mid=1.0 -> sparsity 0.9706,
# mid=3.0 -> sparsity 0.9925. Original mapping had "low"=3.0 (actually more
# severe) and "high"=1.0 (actually less severe). Existing data + param_grid
# relabeled to match; this mapping fixed for consistency on any future
# regeneration. See simulate_splatter.R correction note.
splatter_dropout <- list(
  "none" = list(dropout.type = "none",       dropout.mid = NULL),
  "low"  = list(dropout.type = "experiment", dropout.mid = 1.0),
  "high" = list(dropout.type = "experiment", dropout.mid = 3.0)
)

# NOTE: Sparsity is controlled via bcv.common (biological coefficient of variation),
# NOT mean.rate. Splatter normalizes gene means by their sum during count generation,
# canceling any absolute scaling effect from mean.rate — it has no effect on sparsity.
#
# bcv.common sets negative binomial overdispersion, which directly determines
# structural zero rate. Empirically calibrated at:
#   nGenes = 10000, lib.loc = 7.6 (~2000 UMI/cell), dropout.type = "none"
#   Full table: data/simulated/splatter_calib_bcv2.csv
#
# Achieved structural sparsity at depth = 2000 UMI (lib.loc = 7.6):
#   "0.70" -> bcv.common = 0.01  -> actual = 0.9138
#   "0.80" -> bcv.common = 0.20  -> actual = 0.9179
#   "0.90" -> bcv.common = 0.80  -> actual = 0.9316
#   "0.95" -> bcv.common = 2.00  -> actual = 0.9518
#   "0.98" -> bcv.common = 3.00  -> actual = 0.9623
#
# DESIGN NOTE: Grid labels (0.70 ... 0.98) are ordinal identifiers for five levels
# of structural overdispersion, not absolute sparsity targets. Actual sparsity varies
# jointly with depth and bcv.common. Full range across the depth x bcv.common grid
# spans approximately 0.76 (high depth, low BCV) to 0.99+ (low depth, high BCV).
# Actual achieved sparsity per simulation run is recorded in run metadata at
# simulation time and is the true sparsity variable for all downstream analysis.
splatter_sparsity <- list(
  "0.70" = list(bcv.common = 0.01),
  "0.80" = list(bcv.common = 0.20),
  "0.90" = list(bcv.common = 0.80),
  "0.95" = list(bcv.common = 2.00),
  "0.98" = list(bcv.common = 3.00)
)

# NOTE: lib.loc is log-scale (ln of target depth). lib.scale fixed at 0.2.
# Empirically calibrated at nGenes = 10000, 500 cells, 3 replicates.
# Full table: data/simulated/splatter_calib_depth.csv
splatter_depth <- list(
  "500"   = list(lib.loc = 6.2, lib.scale = 0.2),
  "2000"  = list(lib.loc = 7.6, lib.scale = 0.2),
  "10000" = list(lib.loc = 9.2, lib.scale = 0.2)
)

# -----------------------------------------------------------------------------
# SYMSIM
# -----------------------------------------------------------------------------
# Calibration: calibrate_symsim.R
# Full tables: data/simulated/symsim_calib_sigma.csv
#              data/simulated/symsim_calib_depth.csv
#
# API confirmed:
#   SimulateTrueCounts: randseed=run_id for reproducibility
#   True2ObservedCounts: no seed arg -> set.seed(run_id) before each call
#   No apply_dropout flag; no rangeUMI parameter
#
# SIGMA: Sigma does not meaningfully control sparsity (range 0.62-0.66 across
# Sigma 0.1-2.0 at ngenes=2000). Fixed at 0.4 throughout. Sparsity labels are
# ordinal identifiers, same design as Splatter and scDesign3. Actual achieved
# sparsity recorded per run in metadata.
SYMSIM_SIGMA <- 0.4

# SEPARABILITY: controlled via evf_type and n_de_evf (number of DE EVFs).
# Higher n_de_evf = more genes differ between populations = better separation.
# null uses evf_type="one.population" (single population, no tree structure).
symsim_separability <- list(
  "null"   = list(evf_type = "one.population", n_de_evf = 0L),
  "low"    = list(evf_type = "discrete",        n_de_evf = 2L),
  "medium" = list(evf_type = "discrete",        n_de_evf = 5L),
  "high"   = list(evf_type = "discrete",        n_de_evf = 8L)
)

# DROPOUT: alpha_mean is capture efficiency (higher = more capture = fewer dropouts).
# True2ObservedCounts cannot be disabled; "none" uses alpha_mean=0.97 (near-perfect
# capture, minimal technical zeros beyond count model).
# Confirmed A5 separation: none=0.58 sparsity, low=0.62, high=0.72 at depth_mean=5000.
symsim_dropout <- list(
  "none" = list(alpha_mean = 0.97, alpha_sd = 0.01),
  "low"  = list(alpha_mean = 0.08, alpha_sd = 0.01),
  "high" = list(alpha_mean = 0.02, alpha_sd = 0.01)
)

# DEPTH: depth_mean -> actual UMI/cell is dropout-dependent (sublinear at low
# alpha_mean due to capture saturation). Empirically calibrated per dropout level.
# At high dropout (alpha=0.02), target 10000 UMI is physically unreachable
# (saturates ~4000 UMI at depth_mean=50000); actual UMI recorded in run metadata.
# depth_sd = 10% of depth_mean (min 50).
symsim_depth <- list(
  "none" = list(
    "500"   = list(depth_mean =   500L, depth_sd =   50L),
    "2000"  = list(depth_mean =  2000L, depth_sd =  200L),
    "10000" = list(depth_mean = 10000L, depth_sd = 1000L)
  ),
  "low" = list(
    "500"   = list(depth_mean =   500L, depth_sd =   50L),
    "2000"  = list(depth_mean =  2000L, depth_sd =  200L),
    "10000" = list(depth_mean = 25000L, depth_sd = 2500L)
  ),
  "high" = list(
    "500"   = list(depth_mean =   500L, depth_sd =   50L),
    "2000"  = list(depth_mean =  5000L, depth_sd =  500L),
    "10000" = list(depth_mean = 50000L, depth_sd = 5000L)
  )
)

# -----------------------------------------------------------------------------
# scDESIGN3
# -----------------------------------------------------------------------------

# Separability controlled by cell type subset selected from PBMC 3k reference.
# Reference: TENxPBMCData pbmc3k, top 2000 HVGs, annotated via scran clustering
# and canonical marker genes. Full annotation: data/simulated/pbmc3k_annotated.rds
#
# null   : CD4_T only (1135 cells)           — single population
# low    : CD4_T + CD8_T (1485 cells)        — closely related lymphocytes
# medium : CD4_T + CD8_T + B_cell + NK (2009 cells) — moderate spread
# high   : all 5 types (2695 cells)          — maximally distinct (lymphoid + myeloid)
scdesign3_separability <- list(
  "null"   = list(mu_formula = "1",         cell_types = c("CD4_T")),
  "low"    = list(mu_formula = "cell_type", cell_types = c("CD4_T", "CD8_T")),
  "medium" = list(mu_formula = "cell_type", cell_types = c("CD4_T", "CD8_T", "B_cell", "NK")),
  "high"   = list(mu_formula = "cell_type", cell_types = c("CD4_T", "CD8_T", "B_cell", "NK", "Monocyte"))
)

# pi (zero-inflation) fitted from reference.
# Targeting specific pi requires reference with matching observed zeros.
scdesign3_dropout <- list(
  "none" = list(family_use = "nb",   zero_inflation_pi = NA),
  "low"  = list(family_use = "zinb", zero_inflation_pi = 0.10),
  "high" = list(family_use = "zinb", zero_inflation_pi = 0.40)
)

# Empirically calibrated against PBMC 3k HVG subset (2000 genes).
# Baseline depth of HVG subset = 1116.9 UMI/cell (lower than full matrix
# because counts are distributed across fewer genes after HVG selection).
# Multipliers derived as target / baseline: 500/1117, 2000/1117, 10000/1117.
# Verified: 0.45 -> 500.0, 1.79 -> 1999.9, 8.95 -> 10000.0 UMI/cell.
# Full table: data/simulated/scdesign3_calib_depth.csv
scdesign3_depth <- list(
  "500"   = list(lib_size_multiplier = 0.45),
  "2000"  = list(lib_size_multiplier = 1.79),
  "10000" = list(lib_size_multiplier = 8.95)
)

# =============================================================================
# CALIBRATION FLAGS
# Parameters requiring empirical calibration before simulator runs.
# =============================================================================

calibration_required <- list(
  splatter  = character(0),   # COMPLETE: bcv.common (sparsity) verified in calib_bcv2.csv;
                               #           lib.loc (depth) verified in calib_depth.csv
  symsim    = character(0),   # COMPLETE: Sigma fixed at 0.4 (does not control sparsity at ngenes=2000);
                               #           depth_mean calibrated per dropout level in symsim_calib_depth.csv
  scdesign3 = character(0)   # COMPLETE: depth multipliers verified in scdesign3_calib_depth.csv;
                               #           separability verified via PBMC 3k cluster annotation
)
