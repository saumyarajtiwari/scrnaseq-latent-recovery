# EDA Checkpoint 2.1 — Per-Dataset Technical Summary
# Usage: Rscript technical_summary_one_dataset.R <in_path> <dataset_name> <out_dir>
#
# Computes: overall sparsity (% zeros), mean/median library size (UMIs/cell),
# library size variance, and a Poisson-expected-vs-observed zero-inflation
# index (see script header note below on interpretation).
# Produces: two-panel PNG (library size distribution, genes-detected
# distribution), both log10-scaled given scRNA-seq's typical right skew.
#
# ZERO-INFLATION DEFINITION (explicit, not a universal standard metric):
# For each gene, expected_P(zero) = exp(-mean_count) under a Poisson model;
# observed_P(zero) = actual zero proportion for that gene. Per-dataset
# zero-inflation rate = mean(pmax(observed - expected, 0)) across genes.
# CAVEAT: this is computed on RAW counts (pre-normalization, consistent with
# this project's Step 3 boundary) and will conflate true technical dropout
# with ordinary biological/technical overdispersion, since a Poisson
# comparison has no separate overdispersion parameter. Treat this as a
# relative, cross-dataset comparison metric, not an absolute dropout estimate.

args <- commandArgs(trailingOnly = TRUE)
in_path <- args[1]; dataset_name <- args[2]; out_dir <- args[3]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(in_path)
cts <- assay(sce, "counts")
cat("Loaded", dataset_name, "| dims:", dim(sce), "\n")

lib_size <- Matrix::colSums(cts)
genes_detected <- Matrix::colSums(cts > 0)

overall_sparsity <- 1 - (Matrix::nnzero(cts) / (as.numeric(nrow(cts)) * as.numeric(ncol(cts))))

gene_mean <- Matrix::rowMeans(cts)
gene_obs_zero_rate <- 1 - (Matrix::rowSums(cts > 0) / ncol(cts))
gene_expected_zero_rate <- exp(-gene_mean)
gene_zero_inflation <- pmax(gene_obs_zero_rate - gene_expected_zero_rate, 0)

summary_row <- data.frame(
  dataset_name = dataset_name,
  n_cells = ncol(sce),
  n_genes = nrow(sce),
  overall_sparsity_pct = round(100 * overall_sparsity, 2),
  mean_lib_size = round(mean(lib_size), 2),
  median_lib_size = round(median(lib_size), 2),
  lib_size_variance = round(var(lib_size), 2),
  lib_size_cv = round(sd(lib_size) / mean(lib_size), 4),
  mean_zero_inflation_index = round(mean(gene_zero_inflation), 4),
  median_zero_inflation_index = round(median(gene_zero_inflation), 4)
)

cat("\n=== Summary:", dataset_name, "===\n")
print(summary_row)

summary_path <- file.path(out_dir, paste0(dataset_name, "_technical_summary.csv"))
write.csv(summary_row, summary_path, row.names = FALSE)
cat("Saved summary:", summary_path, "\n")

png(file.path(out_dir, paste0(dataset_name, "_distributions.png")), width = 1000, height = 500)
par(mfrow = c(1, 2))
hist(log10(lib_size + 1), breaks = 50, main = paste0(dataset_name, ": Library Size (log10)"),
     xlab = "log10(UMIs/cell + 1)", col = "steelblue", border = "white")
hist(log10(genes_detected + 1), breaks = 50, main = paste0(dataset_name, ": Genes Detected (log10)"),
     xlab = "log10(genes detected + 1)", col = "darkorange", border = "white")
dev.off()
cat("Saved plot:", file.path(out_dir, paste0(dataset_name, "_distributions.png")), "\n")
