# Fix: scDesign3's dropout "low"/"high" labels were found swapped relative
# to their intended severity, discovered while investigating a non-monotonic
# subspace-recovery-vs-dropout relationship in Step 5.6's plot.
#
# Evidence: achieved_sparsity (actual technical sparsity per file) was
# checked directly against the labeled dropout severity across 6 different
# parameter contexts (varying sparsity, depth, separability). Splatter was
# monotonic (none < low < high) in 6/6 contexts. scDesign3 was non-monotonic
# in 6/6 contexts, with "low" achieved_sparsity consistently EXCEEDING
# "high" achieved_sparsity every single time. Testing the swap hypothesis
# (relabel low<->high) produced perfect monotonicity in all 6/6 contexts,
# confirming a genuine label-swap defect, directly analogous to (though
# independently discovered from) the already-documented Splatter dropout
# label inversion fixed in Step 1.7.
#
# SymSim was also non-monotonic in 5/6 contexts, but the swap hypothesis
# did NOT resolve it there - that issue is left as a documented limitation
# (genuine post-hoc dropout-masking calibration noise), not fixed.
#
# Scope: param_grid.csv itself is NOT modified, since it is shared across
# all three simulators and correctly specifies the intended symmetric
# design (confirmed correct via Splatter's behavior on the same grid rows).
# The fix is applied only to the dropout column for source=="scdesign3"
# rows in the already-joined results tables, where the mislabeling
# originates specifically from how scDesign3's own masking implementation
# applied the "low"/"high" instruction.

swap_scdesign3_dropout <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  idx <- d$source == "scdesign3" & d$dropout %in% c("low", "high")
  cat(path, "- rows swapped:", sum(idx), "\n")
  d$dropout[idx & d$dropout == "low"]  <- "HIGH_TEMP"
  d$dropout[idx & d$dropout == "high"] <- "low"
  d$dropout[idx & d$dropout == "HIGH_TEMP"] <- "high"
  write.csv(d, path, row.names = FALSE)
}

swap_scdesign3_dropout("data/processed/embedding_manifest.csv")
swap_scdesign3_dropout("data/processed/step4_master_results_table.csv")

cat("\nNOTE: step5_axis_sweep_table.csv and all Step 5.2-5.7 plots must be\n")
cat("regenerated after running this fix (dropout is used as a baseline\n")
cat("filter value for the other 5 axis sweeps, not just the dropout sweep\n")
cat("itself), by rerunning scripts/step5_1_organize_by_axis.R and\n")
cat("scripts/step5_2to7_axis_plots.R.\n")
