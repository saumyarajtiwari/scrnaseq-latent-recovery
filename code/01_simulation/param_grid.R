# =============================================================================
# param_grid.R
# Generate the full parameter grid for all simulation runs.
# Output: data/simulated/param_grid.csv
# =============================================================================

library(dplyr, warn.conflicts = FALSE)

# -----------------------------------------------------------------------------
# MAIN GRID
# Full factorial over all technical stress axes.
# n_groups fixed at 5 throughout — not an axis.
# separability = "null" handled exclusively in null control block below.
# -----------------------------------------------------------------------------

main_grid <- expand.grid(
  sparsity        = c(0.70, 0.80, 0.90, 0.95, 0.98),
  depth           = c(500, 2000, 10000),
  dropout         = c("none", "low", "high"),
  separability    = c("low", "medium", "high"),
  n_cells         = c(200, 1000, 5000),
  n_groups        = 5,
  batch           = c("none", "simple", "complex"),
  gene_strategy   = c("all", "hvg2000", "hvg500"),
  clipping        = c("none", "clip99", "log_stabilized"),
  is_null_control = FALSE,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# NULL CONTROL BLOCK
# Single biological population — detects phantom clustering and technical
# artifacts introduced by preprocessing.
# Sparsity varies across all 5 levels; all other axes fixed at baseline.
# Baseline: depth=2000, n_cells=1000, dropout=none, batch=none,
#           gene_strategy=all, clipping=none
# -----------------------------------------------------------------------------

null_grid <- expand.grid(
  sparsity        = c(0.70, 0.80, 0.90, 0.95, 0.98),
  depth           = 2000,
  dropout         = "none",
  separability    = "null",
  n_cells         = 1000,
  n_groups        = 1,
  batch           = "none",
  gene_strategy   = "all",
  clipping        = "none",
  is_null_control = TRUE,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# COMBINE, INDEX, WRITE
# -----------------------------------------------------------------------------

param_grid <- bind_rows(main_grid, null_grid)
param_grid$run_id <- seq_len(nrow(param_grid))
param_grid <- param_grid[, c("run_id", setdiff(names(param_grid), "run_id"))]

cat("=== Parameter Grid Summary ===\n")
cat("Main grid rows   :", nrow(main_grid),  "\n")
cat("Null control rows:", nrow(null_grid),  "\n")
cat("Total rows       :", nrow(param_grid), "\n\n")

cat("Axis levels:\n")
cat("  sparsity     :", paste(sort(unique(main_grid$sparsity)),   collapse = ", "), "\n")
cat("  depth        :", paste(sort(unique(main_grid$depth)),      collapse = ", "), "\n")
cat("  dropout      :", paste(unique(main_grid$dropout),          collapse = ", "), "\n")
cat("  separability :", paste(unique(main_grid$separability),     collapse = ", "), "\n")
cat("  n_cells      :", paste(sort(unique(main_grid$n_cells)),    collapse = ", "), "\n")
cat("  n_groups     : 5 (fixed)\n")
cat("  batch        :", paste(unique(main_grid$batch),            collapse = ", "), "\n")
cat("  gene_strategy:", paste(unique(main_grid$gene_strategy),    collapse = ", "), "\n")
cat("  clipping     :", paste(unique(main_grid$clipping),         collapse = ", "), "\n\n")

out_path <- "data/simulated/param_grid.csv"
write.csv(param_grid, out_path, row.names = FALSE)
cat("Grid written to:", out_path, "\n")
