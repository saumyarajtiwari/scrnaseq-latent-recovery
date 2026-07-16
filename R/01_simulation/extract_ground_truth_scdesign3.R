# =============================================================================
# extract_ground_truth_scdesign3.R
#
# Extracts per-fit_key "true" biological signal (real PBMC3k group-mean
# expression) for all unique scDesign3 fit_keys, using the identical
# reference-subsetting / depth-scaling / subsampling logic as
# simulate_scdesign3.R (same seed = same cells). Since scDesign3 is called
# with return_model=FALSE, the fitted NB-GLM coefficients are not
# accessible; this is a faithful proxy for what that fit targets, computed
# directly from real reference data, not a readout of the fit itself.
#
# separability=="null" rows (single population, e.g. CD4_T only) are
# processed through the same code path as every other row — this
# naturally yields a genes x 1 true_group_means matrix, matching the
# trivial rank-0 ground truth these rows should have.
#
# Read-only w.r.t. param_grid.csv, param_dict.R, pbmc3k_annotated.rds, and
# all of data/simulated/scdesign3/*.rds. Writes only to:
#   data/simulated/ground_truth/scdesign3/scdesign3_truth_run_NNNNN.rds
#     (named by the REPRESENTATIVE run_id of each fit_key group)
#   data/simulated/ground_truth/scdesign3_manifest.csv
# =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

cat("=== scDesign3 Ground-Truth Extraction ===\n")
cat("separability value counts in param_grid.csv:\n")
print(table(param_grid$separability))

fit_cols <- c("depth", "separability", "n_cells", "batch")
fit_key <- function(row) paste(row[fit_cols], collapse="_")
param_grid$fit_key <- apply(param_grid, 1, fit_key)
unique_fits <- unique(param_grid$fit_key)

cat(sprintf("\nUnique fit_keys            : %d\n", length(unique_fits)))
cat(sprintf("null-separability fit_key(s): %s\n\n",
            paste(grep("_null_", unique_fits, value=TRUE), collapse=", ")))

OUT_DIR <- "data/simulated/ground_truth/scdesign3"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("Loading reference...\n")
ref_full <- readRDS("data/simulated/pbmc3k_annotated.rds")
colnames(ref_full) <- paste0("cell_", seq_len(ncol(ref_full)))
counts(ref_full) <- as(counts(ref_full), "dgCMatrix")
cat(sprintf("Reference: %d genes x %d cells\n\n", nrow(ref_full), ncol(ref_full)))

manifest <- data.frame(fit_key=character(), representative_run_id=integer(),
                        n_cell_types=integer(), n_cells_used=integer(),
                        output_path=character(), stringsAsFactors=FALSE)

t_start <- Sys.time()

for (i in seq_along(unique_fits)) {
  fk <- unique_fits[i]
  row1 <- param_grid[param_grid$fit_key == fk, ][1, ]
  run_id <- as.integer(row1$run_id)
  out_path <- file.path(OUT_DIR, sprintf("scdesign3_truth_run_%05d.rds", run_id))

  if (file.exists(out_path)) {
    cat(sprintf("[%2d/%d] SKIP (exists) fit_key=%s\n", i, length(unique_fits), fk))
    manifest <- rbind(manifest, data.frame(fit_key=fk, representative_run_id=run_id,
                       n_cell_types=NA, n_cells_used=NA, output_path=out_path))
    next
  }

  sep_p      <- scdesign3_separability[[row1$separability]]
  cell_types <- sep_p$cell_types
  ref_sub    <- ref_full[, ref_full$cell_type %in% cell_types]

  depth_p    <- scdesign3_depth[[as.character(row1$depth)]]
  cnt_scaled <- round(counts(ref_sub) * depth_p$lib_size_multiplier)

  actual_n <- min(as.integer(row1$n_cells), ncol(ref_sub))
  set.seed(run_id)
  keep_idx <- sample(ncol(ref_sub), actual_n)

  cnt_fit   <- cnt_scaled[, keep_idx]
  types_fit <- as.character(ref_sub$cell_type[keep_idx])
  group_levels <- sort(unique(types_fit))

  true_group_means <- vapply(group_levels, function(g) {
    Matrix::rowMeans(cnt_fit[, types_fit == g, drop=FALSE])
  }, FUN.VALUE = numeric(nrow(cnt_fit)))
  rownames(true_group_means) <- rownames(cnt_fit)
  colnames(true_group_means) <- group_levels

  saveRDS(list(
    true_group_means      = true_group_means,
    fit_key                = fk,
    representative_run_id  = run_id,
    cell_types              = group_levels,
    source                  = "pbmc3k_reference_group_means",
    method_note = paste("Proxy for scDesign3's fitted NB-GLM target mean;",
                         "computed from same subsampled/depth-scaled reference",
                         "cells (return_model=FALSE precludes direct fit access).",
                         "1-column output for separability=='null' rows is the",
                         "expected trivial rank-0 case, not an error.")
  ), out_path, compress=TRUE)

  manifest <- rbind(manifest, data.frame(fit_key=fk, representative_run_id=run_id,
                     n_cell_types=length(group_levels), n_cells_used=actual_n,
                     output_path=out_path))

  cat(sprintf("[%2d/%d] DONE fit_key=%-20s run_id=%05d  types=%d  n=%d\n",
              i, length(unique_fits), fk, run_id, length(group_levels), actual_n))
}

write.csv(manifest, "data/simulated/ground_truth/scdesign3_manifest.csv", row.names=FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
cat(sprintf("\n=== DONE: %d fit_keys in %.1f sec ===\n", length(unique_fits), elapsed))
