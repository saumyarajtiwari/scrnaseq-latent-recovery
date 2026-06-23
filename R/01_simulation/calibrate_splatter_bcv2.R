# =============================================================================
# calibrate_splatter_bcv2.R
# Corrected BCV calibration — nGenes = 10000 to match actual simulation.
#
# Why this supersedes calibrate_splatter_bcv.R:
# At nGenes = 2000 and lib.loc = 7.6 (~2000 UMI/cell), average count per
# gene per cell is 1.0 UMI. At nGenes = 10000, it is 0.2 UMI. This shifts
# the sparsity floor substantially. Calibration must match simulation nGenes.
#
# N_CALIB_CELLS = 200 (reduced for speed — sparsity estimate is stable
# because the count matrix has 200 x 10000 = 2,000,000 entries).
#
# BCV range extended to 3.0 to capture upper targets if achievable.
#
# Output: data/simulated/splatter_calib_bcv2.csv
# =============================================================================

library(splatter)
library(Matrix)

set.seed(42)

N_CALIB_CELLS <- 200
N_CALIB_GENES <- 10000     # matches actual simulation
N_REPS        <- 3

cat("=============================================\n")
cat(" Splatter BCV Calibration v2\n")
cat(" nGenes = 10000 (matches simulation)\n")
cat("=============================================\n\n")
cat(sprintf("Settings: %d cells, %d genes, %d replicates\n\n",
            N_CALIB_CELLS, N_CALIB_GENES, N_REPS))
cat("bcv.common -> % zeros\n")
cat("[lib.loc = 7.6, dropout.type = 'none']\n\n")

bcv_vals    <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00, 3.00)
sparsity_df <- data.frame(bcv_common = bcv_vals, actual_sparsity = NA_real_)

for (i in seq_along(bcv_vals)) {
  bv   <- bcv_vals[i]
  reps <- sapply(seq_len(N_REPS), function(s) {
    sim <- splatSimulate(
      nGenes       = N_CALIB_GENES,
      batchCells   = N_CALIB_CELLS,
      bcv.common   = bv,
      lib.loc      = 7.6,
      lib.scale    = 0.2,
      dropout.type = "none",
      verbose      = FALSE,
      seed         = s * 100
    )
    cnts <- counts(sim)
    sum(cnts == 0L) / length(cnts)
  })
  sparsity_df$actual_sparsity[i] <- round(mean(reps), 4)
  cat(sprintf("  bcv.common = %.2f  ->  sparsity = %.4f\n",
              bv, sparsity_df$actual_sparsity[i]))
}

cat("\nTarget -> nearest bcv.common in tested range:\n")
for (t in c(0.70, 0.80, 0.90, 0.95, 0.98)) {
  idx <- which.min(abs(sparsity_df$actual_sparsity - t))
  cat(sprintf("  target %.2f  ->  bcv.common = %.2f  (actual = %.4f)\n",
              t, sparsity_df$bcv_common[idx], sparsity_df$actual_sparsity[idx]))
}

write.csv(sparsity_df, "data/simulated/splatter_calib_bcv2.csv", row.names = FALSE)
cat("\nTable written to data/simulated/splatter_calib_bcv2.csv\n")
cat("=============================================\n")
