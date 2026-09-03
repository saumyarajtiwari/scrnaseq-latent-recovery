suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

MANIFEST <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)

METHOD_LABELS <- c(
  pca_raw = "Raw PCA",
  pca_libnorm = "Libnorm PCA",
  pca_log = "Log-PCA",
  pca_shiftedlog = "Shifted Log-PCA",
  pca_sctransform_v2 = "Pearson Residual PCA (SCT v2)",
  pca_glmpca = "GLM-PCA"
)
METHOD_ORDER <- names(METHOD_LABELS)
SIMULATORS <- c("scdesign3", "splatter", "symsim")

RUN_BEST  <- 371
RUN_WORST <- 830
RUN_NULL  <- 10940

OUT_BASE <- "figures/eda_checkpoint_3"
for (sub in c("best_case", "worst_case", "null_control")) {
  dir.create(file.path(OUT_BASE, sub), recursive = TRUE, showWarnings = FALSE)
}

get_file <- function(method, run_id, simulator) {
  row <- MANIFEST[MANIFEST$method == method & MANIFEST$run_id == run_id & MANIFEST$source == simulator, ]
  if (nrow(row) != 1) stop("Expected exactly 1 match for ", method, "/", run_id, "/", simulator, ", got ", nrow(row))
  row$file_path[1]
}

make_panel <- function(method, run_id, simulator, title_suffix = "") {
  f <- get_file(method, run_id, simulator)
  d <- readRDS(f)
  df <- data.frame(
    PC1 = d$embedding[, 1],
    PC2 = d$embedding[, 2],
    true_group = as.factor(d$true_group)
  )
  n_groups <- length(unique(df$true_group))
  ggplot(df, aes(x = PC1, y = PC2, color = true_group)) +
    geom_point(size = 0.6, alpha = 0.6) +
    labs(title = paste0(METHOD_LABELS[method], title_suffix), x = "PC1", y = "PC2", color = "Group") +
    theme_minimal(base_size = 9) +
    theme(legend.position = if (n_groups > 1) "right" else "none",
          plot.title = element_text(size = 9, face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
}

cat("=== EDA 3.1: Best-case panels (all 3 simulators) ===\n")
for (sim in SIMULATORS) {
  for (method in METHOD_ORDER) {
    p <- make_panel(method, RUN_BEST, sim, paste0(" - Best-case (", sim, ")"))
    out <- file.path(OUT_BASE, "best_case", paste0(method, "_", sim, "_best.png"))
    ggsave(out, p, width = 5, height = 4, dpi = 150)
  }
  cat("  done:", sim, "\n")
}

cat("\n=== EDA 3.2: Worst-case panels (all 3 simulators) ===\n")
for (sim in SIMULATORS) {
  for (method in METHOD_ORDER) {
    p <- make_panel(method, RUN_WORST, sim, paste0(" - Worst-case (", sim, ")"))
    out <- file.path(OUT_BASE, "worst_case", paste0(method, "_", sim, "_worst.png"))
    ggsave(out, p, width = 5, height = 4, dpi = 150)
  }
  cat("  done:", sim, "\n")
}

cat("\n=== EDA 3.3: Null-control panels (all 3 simulators) ===\n")
for (sim in SIMULATORS) {
  for (method in METHOD_ORDER) {
    p <- make_panel(method, RUN_NULL, sim, paste0(" - Null-control (", sim, ")"))
    out <- file.path(OUT_BASE, "null_control", paste0(method, "_", sim, "_null.png"))
    ggsave(out, p, width = 5, height = 4, dpi = 150)
  }
  cat("  done:", sim, "\n")
}

cat("\n=== EDA 3.4: Combined 6x3 master summary (Splatter) ===\n")
panel_list <- list()
for (method in METHOD_ORDER) {
  panel_list[[length(panel_list) + 1]] <- make_panel(method, RUN_BEST,  "splatter", "\nBest-case")
  panel_list[[length(panel_list) + 1]] <- make_panel(method, RUN_WORST, "splatter", "\nWorst-case")
  panel_list[[length(panel_list) + 1]] <- make_panel(method, RUN_NULL,  "splatter", "\nNull-control")
}
combined <- arrangeGrob(grobs = panel_list, ncol = 3, nrow = 6)
ggsave(file.path(OUT_BASE, "eda3_4_master_summary_splatter.png"), combined,
       width = 15, height = 24, dpi = 150, limitsize = FALSE)

cat("\n=== EDA CHECKPOINT 3 COMPLETE ===\n")
cat("54 individual panels + 1 master summary figure written to", OUT_BASE, "\n")
