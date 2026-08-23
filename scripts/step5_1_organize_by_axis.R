d <- read.csv("data/processed/step4_master_results_table.csv", stringsAsFactors = FALSE)
sim <- d[d$data_type == "simulated", ]

BASELINE <- list(sparsity = 0.9, depth = 2000, dropout = "low", separability = "medium",
                  n_cells = 1000, batch = "none", gene_strategy = "all", clipping = "none")

sweep_axes <- c("sparsity", "depth", "separability", "batch", "dropout", "clipping")

for (axis in sweep_axes) {
  flag_col <- paste0("in_", axis, "_sweep")
  filt <- rep(TRUE, nrow(d))
  for (col in names(BASELINE)) {
    if (col == axis) next
    filt <- filt & (d[[col]] == BASELINE[[col]] | d$data_type == "real")
  }
  d[[flag_col]] <- filt & (d$data_type == "simulated")
}

OUT_PATH <- "data/processed/step5_axis_sweep_table.csv"
write.csv(d, OUT_PATH, row.names = FALSE)

cat("=== STEP 5.1 COMPLETE ===\n")
cat("Baseline used (held fixed except the swept axis):\n")
for (n in names(BASELINE)) cat(" ", n, "=", BASELINE[[n]], "\n")
cat("\nRow counts per sweep:\n")
for (axis in sweep_axes) {
  flag_col <- paste0("in_", axis, "_sweep")
  cat(" ", axis, ":", sum(d[[flag_col]]), "rows\n")
}
