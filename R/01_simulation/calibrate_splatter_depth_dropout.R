# =============================================================================
# calibrate_splatter_depth_dropout.R
#
# CORRECTION (Step 1.7 audit): the original calibrate_splatter.R depth table
# was built at dropout.type="none" only, per its own design note ("Dropout
# layers on top of sparsity") — acknowledging the coupling but never
# quantifying it. This left the Step 1.7 validation comparing dropout=low/
# high files against a dropout-blind depth baseline, producing large
# spurious "deviations" that are actually real, expected dropout-driven
# depth reduction (Splatter's dropout.type="experiment" mechanism zeroes
# real counts, exactly as SymSim's alpha_mean does — already correctly
# calibrated there).
#
# This script extends the same lib.loc sweep to all three dropout levels,
# using the CORRECTED dropout.mid values (low=1.0, high=3.0 — see
# param_dict.R correction note on the low/high label inversion fix).
#
# Output: data/simulated/splatter_calib_depth_dropout.csv
# =============================================================================

library(splatter)
library(Matrix)

set.seed(42)

N_CALIB_CELLS <- 500
N_CALIB_GENES <- 2000
N_REPS        <- 3

lib_loc_vals <- c(5.5, 6.0, 6.2, 6.5, 7.0, 7.6, 8.0, 8.5, 9.0, 9.2)
dropout_configs <- list(
  "none" = list(dropout.type = "none",       dropout.mid = NULL),
  "low"  = list(dropout.type = "experiment", dropout.mid = 1.0),
  "high" = list(dropout.type = "experiment", dropout.mid = 3.0)
)

cat("=============================================\n")
cat(" Splatter Depth x Dropout Calibration\n")
cat("=============================================\n\n")

results <- list()
row_i <- 1L

for (dlabel in names(dropout_configs)) {
  dcfg <- dropout_configs[[dlabel]]
  cat(sprintf("--- dropout = %s ---\n", dlabel))

  for (ll in lib_loc_vals) {
    reps <- sapply(seq_len(N_REPS), function(s) {
      args <- list(
        nGenes = N_CALIB_GENES, batchCells = N_CALIB_CELLS,
        mean.rate = 4.0, lib.loc = ll, lib.scale = 0.2,
        dropout.type = dcfg$dropout.type,
        verbose = FALSE, seed = s * 100
      )
      if (!is.null(dcfg$dropout.mid)) args$dropout.mid <- dcfg$dropout.mid
      sim <- do.call(splatSimulate, args)
      mean(colSums(counts(sim)))
    })
    actual_depth <- round(mean(reps), 1)
    cat(sprintf("  lib.loc = %.1f  ->  mean depth = %8.1f\n", ll, actual_depth))
    results[[row_i]] <- data.frame(dropout_label = dlabel, lib_loc = ll, actual_depth = actual_depth)
    row_i <- row_i + 1L
  }
  cat("\n")
}

depth_dropout_df <- do.call(rbind, results)
write.csv(depth_dropout_df, "data/simulated/splatter_calib_depth_dropout.csv", row.names = FALSE)
cat("Calibration table written to data/simulated/splatter_calib_depth_dropout.csv\n")
