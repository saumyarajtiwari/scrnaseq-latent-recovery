suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

d <- read.csv("data/processed/step5_axis_sweep_table.csv", stringsAsFactors = FALSE)
REP_SIM <- "splatter"

METHOD_LABELS <- c(
  pca_raw = "Raw PCA", pca_libnorm = "Libnorm PCA", pca_log = "Log-PCA",
  pca_shiftedlog = "Shifted Log-PCA", pca_sctransform_v2 = "Pearson Residual PCA (SCT v2)",
  pca_glmpca = "GLM-PCA"
)
METHOD_ORDER <- names(METHOD_LABELS)

BASELINE_OTHER <- list(separability = "medium", n_cells = 1000, batch = "none",
                        gene_strategy = "all", clipping = "none")

get_sweep2d <- function(x_col, y_col, dropout_val = NULL) {
  filt <- d$source == REP_SIM & d$is_null_control == FALSE
  for (col in names(BASELINE_OTHER)) filt <- filt & (d[[col]] == BASELINE_OTHER[[col]])
  if (!is.null(dropout_val)) filt <- filt & (d$dropout == dropout_val)
  d[filt, ]
}

compute_edges <- function(x) {
  # Exact midpoint-based tile edges: each tile spans from the midpoint
  # with its left neighbor to the midpoint with its right neighbor,
  # correctly handling asymmetric spacing (e.g. sparsity 0.9's neighbors
  # at 0.8 and 0.95 are NOT equidistant). Boundary tiles extend by the
  # same gap as their one interior neighbor.
  ux <- sort(unique(x))
  n <- length(ux)
  mids <- (ux[-1] + ux[-n]) / 2
  lo <- c(ux[1] - (mids[1] - ux[1]), mids)
  hi <- c(mids, ux[n] + (ux[n] - mids[n-1]))
  data.frame(value = ux, xmin = lo, xmax = hi)
}

make_heatmap <- function(sub, x_col, y_col, method, title_x, title_y, log_y = FALSE) {
  msub <- sub[sub$method == method, ]
  y_transform <- if (log_y) log10 else identity
  msub$.yval <- y_transform(msub[[y_col]])

  x_edges <- compute_edges(msub[[x_col]])
  y_edges <- compute_edges(msub$.yval)

  msub$xmin <- x_edges$xmin[match(msub[[x_col]], x_edges$value)]
  msub$xmax <- x_edges$xmax[match(msub[[x_col]], x_edges$value)]
  msub$ymin <- y_edges$xmin[match(msub$.yval, y_edges$value)]
  msub$ymax <- y_edges$xmax[match(msub$.yval, y_edges$value)]

  y_breaks <- sort(unique(msub[[y_col]]))
  y_break_positions <- y_transform(y_breaks)

  p <- ggplot(msub, aes(x = .data[[x_col]], y = .yval, z = subspace_recovery_score)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = subspace_recovery_score)) +
    geom_contour(breaks = 0.5, color = "black", linewidth = 0.9) +
    scale_fill_gradient2(low = "darkred", mid = "yellow", high = "darkgreen",
                          midpoint = 0.5, limits = c(0, 1), name = "Recovery\nScore") +
    scale_y_continuous(breaks = y_break_positions, labels = y_breaks) +
    scale_x_continuous(breaks = sort(unique(msub[[x_col]]))) +
    labs(title = METHOD_LABELS[method], x = title_x, y = title_y) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"))
  p
}

cat("=== Sweep 1: Sparsity x Depth (Splatter, 5x3=15 points/method) ===\n")
sweep1 <- get_sweep2d("sparsity", "depth", dropout_val = "low")
cat("Rows:", nrow(sweep1), "(expect", 15*6, ")\n")

panels1 <- list()
for (method in METHOD_ORDER) {
  p <- make_heatmap(sweep1, "sparsity", "depth", method, "Sparsity", "Depth", log_y = TRUE)
  ggsave(file.path("results/step5_phase_space/figures", paste0("sparsity_depth_", method, ".png")),
         p, width = 5, height = 4, dpi = 150)
  panels1[[method]] <- p
}
combined1 <- arrangeGrob(grobs = panels1, ncol = 3, nrow = 2,
                          top = "Step 5.9 Sweep 1: Subspace Recovery - Sparsity x Depth (Splatter)")
ggsave("results/step5_phase_space/figures/sweep1_sparsity_depth_combined.png", combined1,
       width = 15, height = 8, dpi = 150)

cat("\n=== Sweep 2: Dropout x Depth (Splatter, 3x3=9 points/method) ===\n")
filt2 <- d$source == REP_SIM & d$is_null_control == FALSE & d$sparsity == 0.9
for (col in names(BASELINE_OTHER)) filt2 <- filt2 & (d[[col]] == BASELINE_OTHER[[col]])
sweep2 <- d[filt2, ]
sweep2$dropout_num <- match(sweep2$dropout, c("none","low","high"))
cat("Rows:", nrow(sweep2), "(expect", 9*6, ")\n")

panels2 <- list()
for (method in METHOD_ORDER) {
  msub <- sweep2[sweep2$method == method, ]
  msub$.yval <- log10(msub$depth)
  x_edges <- compute_edges(msub$dropout_num)
  y_edges <- compute_edges(msub$.yval)
  msub$xmin <- x_edges$xmin[match(msub$dropout_num, x_edges$value)]
  msub$xmax <- x_edges$xmax[match(msub$dropout_num, x_edges$value)]
  msub$ymin <- y_edges$xmin[match(msub$.yval, y_edges$value)]
  msub$ymax <- y_edges$xmax[match(msub$.yval, y_edges$value)]
  y_breaks <- sort(unique(msub$depth))
  p <- ggplot(msub, aes(x = dropout_num, y = .yval, z = subspace_recovery_score)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = subspace_recovery_score)) +
    geom_contour(breaks = 0.5, color = "black", linewidth = 0.9) +
    scale_fill_gradient2(low = "darkred", mid = "yellow", high = "darkgreen",
                          midpoint = 0.5, limits = c(0, 1), name = "Recovery\nScore") +
    scale_x_continuous(breaks = 1:3, labels = c("none","low","high")) +
    scale_y_continuous(breaks = log10(y_breaks), labels = y_breaks) +
    labs(title = METHOD_LABELS[method], x = "Dropout", y = "Depth") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"))
  ggsave(file.path("results/step5_phase_space/figures", paste0("dropout_depth_", method, ".png")),
         p, width = 5, height = 4, dpi = 150)
  panels2[[method]] <- p
}
combined2 <- arrangeGrob(grobs = panels2, ncol = 3, nrow = 2,
                          top = "Step 5.9 Sweep 2: Subspace Recovery - Dropout x Depth (Splatter)")
ggsave("results/step5_phase_space/figures/sweep2_dropout_depth_combined.png", combined2,
       width = 15, height = 8, dpi = 150)

cat("\n=== STEP 5.9 COMPLETE ===\n")
