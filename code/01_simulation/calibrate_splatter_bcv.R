# =============================================================================
# calibrate_splatter_bcv.R
# Calibrate bcv.common -> actual % zeros at dropout.type = "none".
# bcv.common (biological coefficient of variation) controls NB overdispersion,
# which directly sets structural sparsity.
# mean.rate does not control sparsity (Splatter normalizes gene means by their
# sum, so absolute scale cancels). This script replaces that approach.
# Output: data/simulated/splatter_calib_bcv.csv
# =============================================================================

library(splatter)
library(Matrix)

set.seed(42)

N_CALIB_CELLS <- 500
N_CALIB_GENES <- 2000
N_REPS        <- 3

cat("=============================================\n")
cat(" Splatter BCV -> Sparsity Calibration\n")
cat("=============================================\n\n")
cat(sprintf("Settings: %d cells, %d genes, %d replicates\n\n",
            N_CALIB_CELLS, N_CALIB_GENES, N_REPS))
cat("PART 1 — bcv.common -> % zeros\n")
cat("     [lib.loc = 7.6, dropout.type = 'none' fixed]\n\n")

bcv_vals    <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00)
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

cat("\n  Target -> recommended bcv.common:\n")
for (t in c(0.70, 0.80, 0.90, 0.95, 0.98)) {
  idx <- which.min(abs(sparsity_df$actual_sparsity - t))
  cat(sprintf("  target %.2f  ->  bcv.common = %.2f  (actual = %.4f)\n",
              t, sparsity_df$bcv_common[idx], sparsity_df$actual_sparsity[idx]))
}

write.csv(sparsity_df, "data/simulated/splatter_calib_bcv.csv", row.names = FALSE)
cat("\n  Table written to data/simulated/splatter_calib_bcv.csv\n")
cat("=============================================\n")
