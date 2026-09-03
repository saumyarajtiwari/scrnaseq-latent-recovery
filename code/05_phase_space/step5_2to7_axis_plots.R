suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

d <- read.csv("data/processed/step5_axis_sweep_table.csv", stringsAsFactors = FALSE)

METHOD_LABELS <- c(
  pca_raw = "Raw PCA", pca_libnorm = "Libnorm PCA", pca_log = "Log-PCA",
  pca_shiftedlog = "Shifted Log-PCA", pca_sctransform_v2 = "Pearson Residual PCA (SCT v2)",
  pca_glmpca = "GLM-PCA"
)
d$method_label <- METHOD_LABELS[d$method]

AXIS_SPECS <- list(
  sparsity     = list(flag = "in_sparsity_sweep",     x = "sparsity",     order = c(0.7,0.8,0.9,0.95,0.98), title = "Sparsity", step = "5.2"),
  depth        = list(flag = "in_depth_sweep",        x = "depth",        order = c(500,2000,10000),        title = "Sequencing Depth", step = "5.3"),
  separability = list(flag = "in_separability_sweep", x = "separability", order = c("low","medium","high"), title = "Cell-Type Separability", step = "5.4"),
  batch        = list(flag = "in_batch_sweep",        x = "batch",        order = c("none","simple","complex"), title = "Batch Complexity", step = "5.5"),
  dropout      = list(flag = "in_dropout_sweep",      x = "dropout",      order = c("none","low","high"),   title = "Dropout Rate", step = "5.6"),
  clipping     = list(flag = "in_clipping_sweep",     x = "clipping",     order = c("none","log_stabilized","clip99"), title = "Clipping Scale", step = "5.7")
)

for (axis_name in names(AXIS_SPECS)) {
  spec <- AXIS_SPECS[[axis_name]]
  sub <- d[d[[spec$flag]] == TRUE, ]
  sub[[spec$x]] <- factor(sub[[spec$x]], levels = spec$order)

  p <- ggplot(sub, aes(x = .data[[spec$x]], y = subspace_recovery_score,
                        color = method_label, group = method_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    facet_wrap(~ source, ncol = 3) +
    labs(title = paste0("Step ", spec$step, ": Subspace Recovery vs. ", spec$title),
         subtitle = "Dashed line = 0.5 failure threshold. Faceted by simulator.",
         x = spec$title, y = "Subspace Recovery Score", color = "Method") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  out_file <- file.path("figures/step5_axis_sweeps", paste0("step", gsub("\\.", "_", spec$step), "_", axis_name, ".png"))
  ggsave(out_file, p, width = 12, height = 4.5, dpi = 150)
  cat("Saved:", out_file, "\n")
}

cat("\n=== STEPS 5.2-5.7 COMPLETE ===\n")
