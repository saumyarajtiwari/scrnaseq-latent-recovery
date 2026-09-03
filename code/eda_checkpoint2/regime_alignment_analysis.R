# EDA 2.4 — Real vs. Simulated Regime Alignment: nearest-neighbor match in
# (log10 depth, achieved sparsity %) space, per real dataset, against all
# non-null-control simulated runs. For datasets outside the simulated depth
# range (per EDA 2.2), the "nearest" match is necessarily an edge/extrapolated
# match, not a true interpolated one -- flagged explicitly, not glossed over.

sim <- read.csv("results/eda_checkpoint2/simulated_phase_space_summary.csv", stringsAsFactors = FALSE)
real <- read.csv("results/eda_checkpoint2/all_datasets_technical_summary.csv", stringsAsFactors = FALSE)
sim_non_null <- sim[!sim$is_null_control, ]

sim_log_depth <- log10(sim_non_null$mean_lib_size)
sim_sparsity <- sim_non_null$achieved_sparsity_pct

# Normalize both axes to comparable scale before computing distance
depth_range <- range(c(sim_log_depth, log10(real$mean_lib_size)))
sparsity_range <- range(c(sim_sparsity, real$overall_sparsity_pct))

norm_depth_sim <- (sim_log_depth - depth_range[1]) / diff(depth_range)
norm_sparsity_sim <- (sim_sparsity - sparsity_range[1]) / diff(sparsity_range)

results <- list()
for (i in seq_len(nrow(real))) {
  d <- real[i, ]
  norm_depth_real <- (log10(d$mean_lib_size) - depth_range[1]) / diff(depth_range)
  norm_sparsity_real <- (d$overall_sparsity_pct - sparsity_range[1]) / diff(sparsity_range)

  dist <- sqrt((norm_depth_sim - norm_depth_real)^2 + (norm_sparsity_sim - norm_sparsity_real)^2)
  nearest_idx <- order(dist)[1:5]
  nearest <- sim_non_null[nearest_idx, ]

  is_extrapolated <- log10(d$mean_lib_size) > max(sim_log_depth) | log10(d$mean_lib_size) < min(sim_log_depth)

  cat("\n===", d$dataset_name, "===\n")
  cat("Real: log10(depth)=", round(log10(d$mean_lib_size),2), "sparsity%=", d$overall_sparsity_pct, "\n")
  cat("Extrapolated (outside simulated depth range):", is_extrapolated, "\n")
  cat("Top 5 nearest simulated regimes (by simulator, sparsity_label, depth_label, dropout, separability, distance):\n")
  print(data.frame(nearest[, c("simulator","sparsity_label","depth_label","dropout","separability")],
                    distance = round(dist[nearest_idx], 4)))

  results[[d$dataset_name]] <- list(real = d, nearest = nearest, dist = dist[nearest_idx], is_extrapolated = is_extrapolated)
}

saveRDS(results, "results/eda_checkpoint2/regime_alignment_results.rds")
cat("\nSaved: results/eda_checkpoint2/regime_alignment_results.rds\n")
