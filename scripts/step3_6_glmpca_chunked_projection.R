suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(glmpca)
})

MAX_GENES     <- 3000
PCA_NV        <- 30
SEED          <- 42
TRAIN_SIZE    <- 2000
TRAIN_LR      <- 0.0008
OUT_ROOT      <- "data/processed/pca_glmpca/real"

# Baron validation reference (documented, not re-derived here):
# ARI projected=0.421 vs full-fit=0.652; Silhouette projected=0.121 vs full-fit=0.073
# (Baron, held-out 6569 cells, standardized embeddings, k-means k=14)
VALIDATION_NOTE <- paste(
  "Fitted via chunked fixed-loadings projection: GLM-PCA trained on a",
  TRAIN_SIZE, "-cell representative subsample, remaining cells projected",
  "onto fixed loadings/dispersion via per-cell deviance minimization",
  "(BFGS, using glmpca's own dev_func/infograd). Validated against Baron",
  "(the only real dataset with both a native full-data GLM-PCA fit and a",
  "feasible held-out test): standardized-embedding downstream metrics on",
  "6569 held-out cells showed ARI=0.421 (projected) vs 0.652 (native full",
  "fit), both evaluated with the project's own k-means/ARI protocol. This",
  "is a real, quantified degradation vs. a native full fit -- reported",
  "transparently since no native full fit is computable for this file",
  "(65,690/61,292 cells exceeds available memory; see Step 3.6 GLM-PCA",
  "documentation). Use with this caveat in downstream analysis."
)

project_one_cell <- function(y_vec, lib_size, V, coefX_vec, offset_mean, gf, L) {
  offset_i <- log(lib_size) - offset_mean
  deviance_fn <- function(u) gf$dev_func(y_vec, offset_i + coefX_vec + as.numeric(V %*% u))
  grad_fn <- function(u) {
    eta <- offset_i + coefX_vec + as.numeric(V %*% u)
    ig <- gf$infograd(y_vec, eta)
    as.numeric(-2 * t(V) %*% ig$grad)
  }
  optim(par = rep(0, L), fn = deviance_fn, gr = grad_fn, method = "BFGS", control = list(maxit = 100))$par
}

run_chunked_projection <- function(sce_path, dataset_name) {
  cat("=== ", dataset_name, " ===\n")
  t0 <- Sys.time()
  sce <- readRDS(sce_path)
  counts_mat <- counts(sce)
  if (is.null(rownames(counts_mat))) rownames(counts_mat) <- paste0("gene", seq_len(nrow(counts_mat)))
  if (is.null(colnames(counts_mat))) colnames(counts_mat) <- paste0("cell", seq_len(ncol(counts_mat)))

  gene_totals <- Matrix::rowSums(counts_mat)
  counts_mat <- counts_mat[gene_totals > 0, , drop = FALSE]

  n_genes_orig <- nrow(counts_mat)
  hvg_capped <- FALSE
  if (nrow(counts_mat) > MAX_GENES) {
    gene_means <- Matrix::rowMeans(counts_mat)
    top_genes <- order(gene_means, decreasing = TRUE)[1:MAX_GENES]
    counts_mat <- counts_mat[top_genes, , drop = FALSE]
    hvg_capped <- TRUE
  }

  n_cells <- ncol(counts_mat)
  set.seed(SEED)
  train_idx <- sample(seq_len(n_cells), min(TRAIN_SIZE, floor(n_cells * 0.5)))
  test_idx <- setdiff(seq_len(n_cells), train_idx)
  cat("  Training on", length(train_idx), "cells, projecting", length(test_idx), "cells\n")

  train_mat <- as.matrix(counts_mat[, train_idx])
  set.seed(SEED)
  fit <- glmpca(train_mat, L = PCA_NV, fam = "nb",
                ctl = list(verbose = FALSE, maxIter = 25, minIter = 2, lr = TRAIN_LR))

  V <- as.matrix(fit$loadings)
  coefX_vec <- as.numeric(as.matrix(fit$coefX)[, 1])
  offset_mean <- mean(log(Matrix::colSums(train_mat)))
  gf <- fit$glmpca_family

  test_mat <- as.matrix(counts_mat[, test_idx])
  test_lib_sizes <- Matrix::colSums(test_mat)
  cat("  Projecting", ncol(test_mat), "cells...\n")
  projected_U <- t(sapply(seq_len(ncol(test_mat)), function(i) {
    project_one_cell(test_mat[, i], test_lib_sizes[i], V, coefX_vec, offset_mean, gf, PCA_NV)
  }))
  rownames(projected_U) <- colnames(test_mat)

  train_factors <- as.matrix(fit$factors)
  rownames(train_factors) <- colnames(train_mat)

  full_embedding <- rbind(train_factors, projected_U)
  full_embedding <- full_embedding[colnames(counts_mat), ]  # restore original cell order

  out <- list(
    embedding = full_embedding,
    loadings = V,
    nv_used = PCA_NV,
    seed = SEED,
    method = "glmpca_chunked_projection",
    family = "negative_binomial",
    link = "log",
    n_genes_original = n_genes_orig,
    n_genes_used = nrow(counts_mat),
    hvg_capped = hvg_capped,
    train_size = length(train_idx),
    projected_size = length(test_idx),
    validation_note = VALIDATION_NOTE,
    source_file = sce_path,
    true_group = sce$true_group,
    batch_id = sce$batch_id,
    runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  out_dir <- file.path(OUT_ROOT, dataset_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(basename(tools::file_path_sans_ext(sce_path)), "_glmpca.rds"))
  tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
  saveRDS(out, tmp_file)
  file.rename(tmp_file, out_file)

  cat("  Done. Total runtime:", round(out$runtime_sec, 1), "sec. Saved to:", out_file, "\n\n")
  out_file
}

run_chunked_projection("data/real/pbmc68k/pbmc68k_harmonized.rds", "pbmc68k")
run_chunked_projection("data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds", "tabula_sapiens_lung")

cat("=== CHUNKED PROJECTION COMPLETE ===\n")
