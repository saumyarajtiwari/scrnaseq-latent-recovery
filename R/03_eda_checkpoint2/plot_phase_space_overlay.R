# EDA Checkpoint 2.2 (part 2) — Phase-space overlay: real datasets vs.
# simulated parameter grid, sparsity (%) vs. depth (mean UMIs/cell, log10).

sim <- read.csv("results/eda_checkpoint2/simulated_phase_space_summary.csv", stringsAsFactors = FALSE)
real <- read.csv("results/eda_checkpoint2/all_datasets_technical_summary.csv", stringsAsFactors = FALSE)

cat("Simulated rows:", nrow(sim), "| Real rows:", nrow(real), "\n")

sim_non_null <- sim[!sim$is_null_control, ]

png("results/eda_checkpoint2/phase_space_overlay.png", width = 1100, height = 850)

plot(log10(sim_non_null$mean_lib_size), sim_non_null$achieved_sparsity_pct,
     pch = 16, col = adjustcolor("gray60", alpha.f = 0.15), cex = 0.6,
     xlab = "log10(mean UMIs / cell)", ylab = "Achieved sparsity (%)",
     main = "Phase-Space Overlay: Real Datasets vs. Simulated Parameter Grid",
     xlim = range(c(log10(sim_non_null$mean_lib_size), log10(real$mean_lib_size))),
     ylim = c(0, 100))

points(log10(real$mean_lib_size), real$overall_sparsity_pct,
       pch = 17, col = "red3", cex = 2.2)

text(log10(real$mean_lib_size), real$overall_sparsity_pct,
     labels = real$dataset_name, pos = 3, col = "red3", font = 2, cex = 0.9)

legend("bottomleft",
       legend = c("Simulated runs (all simulators, non-null)", "Real datasets"),
       col = c("gray60", "red3"), pch = c(16, 17), pt.cex = c(1, 1.5), bty = "n")

dev.off()
cat("Saved: results/eda_checkpoint2/phase_space_overlay.png\n")

# Coverage check: is each real dataset within the simulated grid's convex range?
cat("\n=== Simulated grid range ===\n")
cat("log10(mean_lib_size): [", round(min(log10(sim_non_null$mean_lib_size)),2), ",",
    round(max(log10(sim_non_null$mean_lib_size)),2), "]\n")
cat("achieved_sparsity_pct: [", round(min(sim_non_null$achieved_sparsity_pct),2), ",",
    round(max(sim_non_null$achieved_sparsity_pct),2), "]\n\n")

for (i in seq_len(nrow(real))) {
  d <- real[i, ]
  log_depth <- log10(d$mean_lib_size)
  in_depth_range <- log_depth >= min(log10(sim_non_null$mean_lib_size)) &&
                    log_depth <= max(log10(sim_non_null$mean_lib_size))
  in_sparsity_range <- d$overall_sparsity_pct >= min(sim_non_null$achieved_sparsity_pct) &&
                       d$overall_sparsity_pct <= max(sim_non_null$achieved_sparsity_pct)
  cat(d$dataset_name, "- depth in range:", in_depth_range, "| sparsity in range:", in_sparsity_range, "\n")
}
