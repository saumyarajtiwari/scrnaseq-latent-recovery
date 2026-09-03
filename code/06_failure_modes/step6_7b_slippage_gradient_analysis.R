# Step 6.7b — Formal monotonicity check and first-crossing-threshold
# identification, per project convention: monotonicity checked on MEANS
# aggregated across replicates, not individual replicate values.

d <- read.csv("data/processed/step6_7_subspace_rotation_slippage.csv", stringsAsFactors = FALSE)

# Ordinal stressor axes with their theoretically expected direction of
# increasing stress (used only for reporting orientation, not to force a result)
ordinal_axes <- list(
  sparsity = list(levels = c("0.7","0.8","0.9","0.95","0.98"), col = "sparsity"),
  depth = list(levels = c("500","2000","10000"), col = "depth"),  # more depth = LESS stress
  dropout = list(levels = c("none","low","high"), col = "dropout"),
  separability = list(levels = c("high","medium","low"), col = "separability"),  # lower sep = MORE stress
  batch = list(levels = c("none","simple","complex"), col = "batch")
)

methods <- unique(d$method)
results <- list()

for (axis_name in names(ordinal_axes)) {
  axis <- ordinal_axes[[axis_name]]
  for (m in methods) {
    for (pc in c("angle_pc1", "angle_pc2", "angle_pc3")) {
      sub <- d[d$method == m, c(axis$col, pc)]
      colnames(sub) <- c("level", "angle")
      sub$level <- factor(as.character(sub$level), levels = axis$levels, ordered = TRUE)
      sub <- sub[!is.na(sub$level) & !is.na(sub$angle), ]
      if (nrow(sub) == 0) next

      agg <- aggregate(angle ~ level, data = sub, FUN = mean)
      agg <- agg[order(agg$level), ]
      if (nrow(agg) < 2) next

      level_numeric <- as.numeric(agg$level)
      rho <- suppressWarnings(cor(level_numeric, agg$angle, method = "spearman"))
      is_monotonic_increasing <- all(diff(agg$angle) >= -1e-6)
      is_monotonic_decreasing <- all(diff(agg$angle) <= 1e-6)

      first_cross <- agg$level[which(agg$angle > 30)]
      first_cross_level <- if (length(first_cross) > 0) as.character(first_cross[1]) else "never_crosses_30"
      always_above_30 <- all(agg$angle > 30)

      results[[length(results) + 1]] <- data.frame(
        axis = axis_name, method = m, pc = pc, spearman_rho = rho,
        monotonic_increasing = is_monotonic_increasing,
        monotonic_decreasing = is_monotonic_decreasing,
        first_level_crossing_30deg = first_cross_level,
        always_above_30deg_ceiling_effect = always_above_30,
        min_angle = min(agg$angle), max_angle = max(agg$angle)
      )
    }
  }
}

results_df <- do.call(rbind, results)
write.csv(results_df, "data/processed/step6_7b_slippage_gradient_analysis.csv", row.names = FALSE)

cat("=== Cleanest monotonic axes (|rho| >= 0.9, PC1 only) ===\n")
pc1_results <- results_df[results_df$pc == "angle_pc1", ]
print(pc1_results[abs(pc1_results$spearman_rho) >= 0.9,
                   c("axis","method","spearman_rho","monotonic_increasing","monotonic_decreasing")])

cat("\n=== Ceiling-effect cases (always above 30deg regardless of stressor level) ===\n")
print(unique(pc1_results[pc1_results$always_above_30deg_ceiling_effect,
                          c("axis","method")]))

cat("\n=== Non-monotonic cases (PC1, worth flagging) ===\n")
print(pc1_results[!pc1_results$monotonic_increasing & !pc1_results$monotonic_decreasing,
                   c("axis","method","spearman_rho","min_angle","max_angle")])

cat("\nWritten to data/processed/step6_7b_slippage_gradient_analysis.csv\n")
