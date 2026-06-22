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

splatter_dropout <- list(
  "none" = list(dropout.type = "none",       dropout.mid = NULL),
  "low"  = list(dropout.type = "experiment", dropout.mid = 3.0),
  "high" = list(dropout.type = "experiment", dropout.mid = 1.0)
)

# NOTE: mean.rate values are starting estimates.
# Empirical calibration required before full grid execution.
splatter_sparsity <- list(
  "0.70" = list(mean.rate = 7.0),
  "0.80" = list(mean.rate = 4.0),
  "0.90" = list(mean.rate = 2.0),
  "0.95" = list(mean.rate = 1.0),
  "0.98" = list(mean.rate = 0.3)
)

# NOTE: lib.loc is log-scale (ln of target depth). lib.scale fixed at 0.2.
# Starting estimates — empirical calibration required.
splatter_depth <- list(
  "500"   = list(lib.loc = 6.2, lib.scale = 0.2),
  "2000"  = list(lib.loc = 7.6, lib.scale = 0.2),
  "10000" = list(lib.loc = 9.2, lib.scale = 0.2)
)

# -----------------------------------------------------------------------------
# SYMSIM
# -----------------------------------------------------------------------------

symsim_separability <- list(
  "null"   = list(n_cell_types = 1, evf_center_sd = 0.0),
  "low"    = list(n_cell_types = 5, evf_center_sd = 0.5),
  "medium" = list(n_cell_types = 5, evf_center_sd = 1.0),
  "high"   = list(n_cell_types = 5, evf_center_sd = 2.0)
)

# NOTE: alpha_mean is capture efficiency.
# Higher alpha_mean = more capture = fewer dropouts.
# low dropout  → alpha_mean = 0.08 (high capture)
# high dropout → alpha_mean = 0.02 (low capture)
symsim_dropout <- list(
  "none" = list(apply_dropout = FALSE, alpha_mean = NULL, alpha_sd = NULL),
  "low"  = list(apply_dropout = TRUE,  alpha_mean = 0.08, alpha_sd = 0.02),
  "high" = list(apply_dropout = TRUE,  alpha_mean = 0.02, alpha_sd = 0.01)
)

# NOTE: rangeUMI brackets require empirical verification.
symsim_depth <- list(
  "500"   = list(rangeUMI = c(200,  800)),
  "2000"  = list(rangeUMI = c(1000, 3000)),
  "10000" = list(rangeUMI = c(6000, 14000))
)

# NOTE: Sigma/bimod interaction for sparsity targeting requires
# empirical calibration — same calibration pass as Splatter mean.rate.

# -----------------------------------------------------------------------------
# scDESIGN3
# -----------------------------------------------------------------------------

# Separability controlled by mu_formula + reference dataset.
# Reference datasets specified in scDesign3 simulation script.
scdesign3_separability <- list(
  "null"   = list(mu_formula = "1",         reference_structure = "single_population"),
  "low"    = list(mu_formula = "cell_type", reference_structure = "overlapping_clusters"),
  "medium" = list(mu_formula = "cell_type", reference_structure = "moderate_separation"),
  "high"   = list(mu_formula = "cell_type", reference_structure = "well_separated")
)

# pi (zero-inflation) fitted from reference.
# Targeting specific pi requires reference with matching observed zeros.
scdesign3_dropout <- list(
  "none" = list(family_use = "nb",   zero_inflation_pi = NA),
  "low"  = list(family_use = "zinb", zero_inflation_pi = 0.10),
  "high" = list(family_use = "zinb", zero_inflation_pi = 0.40)
)

# Assumes reference baseline depth ~2000 UMI.
# Multiplier verified at scDesign3 script stage.
scdesign3_depth <- list(
  "500"   = list(lib_size_multiplier = 0.25),
  "2000"  = list(lib_size_multiplier = 1.00),
  "10000" = list(lib_size_multiplier = 5.00)
)

# =============================================================================
# CALIBRATION FLAGS
# Parameters below are starting estimates requiring empirical calibration
# before full grid execution.
# =============================================================================

calibration_required <- list(
  splatter  = c("mean.rate (sparsity)", "lib.loc (depth)"),
  symsim    = c("Sigma/bimod (sparsity)", "rangeUMI (depth)"),
  scdesign3 = c("lib_size_multiplier (depth)", "zero_inflation_pi via reference selection")
)
