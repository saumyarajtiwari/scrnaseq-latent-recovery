# =============================================================================
# calibrate_splatter.R
# Empirically map Splatter parameters to target sparsity and depth values.
# Must be run before simulate_splatter.R.
#
# Design note: Sparsity is calibrated at dropout.type = "none".
# This isolates structural zeros (count model) from technical zeros (dropout),
# keeping the two axes independent. Dropout layers on top of sparsity.
#
# Output: data/simulated/splatter_calib_sparsity.csv
#         data/simulated/splatter_calib_depth.csv
# =============================================================================

library(splatter)
library(Matrix)

set.seed(42)

N_CALIB_CELLS <- 500
N_CALIB_GENES <- 2000
N_REPS        <- 3

cat("=============================================\n")
cat(" Splatter Calibration Run\n")
cat("=============================================\n\n")
cat(sprintf("Settings: %d cells, %d genes, %d replicates per value\n\n",
            N_CALIB_CELLS, N_CALIB_GENES, N_REPS))

# -----------------------------------------------------------------------------
# PART 1 — Sparsity: mean.rate → actual % zeros
# lib.loc fixed at 7.6 (baseline ~2000 UMI depth)
# dropout.type = "none" — structural zeros only
# -----------------------------------------------------------------------------

cat("PART 1 — Sparsity calibration (mean.rate -> % zeros)\n")
cat("     [dropout.type = 'none', lib.loc = 7.6 fixed]\n\n")

mean_rate_vals <- c(0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 7.0, 10.0, 15.0, 20.0)
sparsity_df    <- data.frame(mean_rate = mean_rate_vals, actual_sparsity = NA_real_)

for (i in seq_along(mean_rate_vals)) {
  mr   <- mean_rate_vals[i]
  reps <- sapply(seq_len(N_REPS), function(s) {
    sim <- splatSimulate(
      nGenes       = N_CALIB_GENES,
      batchCells   = N_CALIB_CELLS,
      mean.rate    = mr,
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
  cat(sprintf("  mean.rate = %5.1f  ->  sparsity = %.4f\n",
              mr, sparsity_df$actual_sparsity[i]))
}

cat("\n  Target -> recommended mean.rate:\n")
for (t in c(0.70, 0.80, 0.90, 0.95, 0.98)) {
  idx <- which.min(abs(sparsity_df$actual_sparsity - t))
  cat(sprintf("  target %.2f  ->  mean.rate = %5.1f  (actual = %.4f)\n",
              t, sparsity_df$mean_rate[idx], sparsity_df$actual_sparsity[idx]))
}

# -----------------------------------------------------------------------------
# PART 2 — Depth: lib.loc → mean UMIs per cell
# mean.rate fixed at 4.0 (moderate baseline)
# -----------------------------------------------------------------------------

cat("\nPART 2 — Depth calibration (lib.loc -> mean UMIs/cell)\n")
cat("     [mean.rate = 4.0, dropout.type = 'none' fixed]\n\n")

lib_loc_vals <- c(5.5, 6.0, 6.2, 6.5, 7.0, 7.6, 8.0, 8.5, 9.0, 9.2)
depth_df     <- data.frame(lib_loc = lib_loc_vals, actual_depth = NA_real_)

for (i in seq_along(lib_loc_vals)) {
  ll   <- lib_loc_vals[i]
  reps <- sapply(seq_len(N_REPS), function(s) {
    sim <- splatSimulate(
      nGenes       = N_CALIB_GENES,
      batchCells   = N_CALIB_CELLS,
      mean.rate    = 4.0,
      lib.loc      = ll,
      lib.scale    = 0.2,
      dropout.type = "none",
      verbose      = FALSE,
      seed         = s * 100
    )
    mean(colSums(counts(sim)))
  })
  depth_df$actual_depth[i] <- round(mean(reps), 1)
  cat(sprintf("  lib.loc = %.1f  ->  mean depth = %8.1f\n", ll, depth_df$actual_depth[i]))
}

cat("\n  Target -> recommended lib.loc:\n")
for (d in c(500, 2000, 10000)) {
  idx <- which.min(abs(depth_df$actual_depth - d))
  cat(sprintf("  target %5d  ->  lib.loc = %.1f  (actual = %.1f)\n",
              d, depth_df$lib_loc[idx], depth_df$actual_depth[idx]))
}

# -----------------------------------------------------------------------------
# WRITE OUTPUT
# -----------------------------------------------------------------------------

write.csv(sparsity_df, "data/simulated/splatter_calib_sparsity.csv", row.names = FALSE)
write.csv(depth_df,    "data/simulated/splatter_calib_depth.csv",    row.names = FALSE)
cat("\n  Calibration tables written to data/simulated/\n")
cat("=============================================\n")
cat(" Use the recommended values above to update\n")
cat(" splatter_sparsity and splatter_depth in\n")
cat(" R/01_simulation/param_dict.R\n")
cat("=============================================\n")
