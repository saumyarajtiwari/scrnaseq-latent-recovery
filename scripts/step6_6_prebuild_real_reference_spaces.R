# Sequentially pre-build the reference space cache for all 6 real datasets,
# single-process (no mclapply), with explicit intermediate cleanup between
# steps. Found necessary after 4/500 dry-run worker crashes concentrated on
# the two highest-gene-count real datasets (Tabula Sapiens Lung: 56,139
# genes; Tasic2018: 42,865 genes) -- multiple full-size sparse matrix copies
# (normalized, logged, logged^2) held simultaneously per worker, times up to
# 4 concurrent workers on different large files, plausibly exceeded
# available memory. This pre-build removes that risk entirely: once cached,
# the main parallel run only ever reads the cache for real data.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(BiocNeighbors)
  library(irlba)
})

N_HVG <- 500
N_PCS_REF <- 15
K_NEIGHBORS <- 15

ref_cache_dir <- "data/processed/step6_reference_space_cache"
dir.create(ref_cache_dir, showWarnings = FALSE, recursive = TRUE)

inv <- read.csv("data/real_data_inventory.csv", stringsAsFactors = FALSE)

for (i in seq_len(nrow(inv))) {
  source_file <- inv$file_path[i]
  cache_key <- gsub("[/\\\\]", "_", source_file)
  cache_path <- file.path(ref_cache_dir, paste0(cache_key, ".rds"))

  if (file.exists(cache_path)) {
    cat(sprintf("[%s] already cached, skipping\n", inv$dataset_name[i]))
    next
  }

  cat(sprintf("[%s] building reference space...\n", inv$dataset_name[i]))
  t0 <- Sys.time()

  src <- readRDS(source_file)
  counts_mat <- counts(src)
  rm(src); gc(verbose = FALSE)

  lib_sizes <- Matrix::colSums(counts_mat)
  safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
  normalized <- counts_mat %*% Diagonal(x = 10000 / safe_lib)
  rm(counts_mat); gc(verbose = FALSE)

  logged <- log1p(normalized)
  rm(normalized); gc(verbose = FALSE)

  n_cells_local <- ncol(logged)
  gene_mean <- Matrix::rowMeans(logged)
  logged_sq <- logged^2
  gene_mean_sq <- Matrix::rowMeans(logged_sq)
  rm(logged_sq); gc(verbose = FALSE)
  gene_var <- (gene_mean_sq - gene_mean^2) * (n_cells_local / (n_cells_local - 1))

  n_hvg <- min(N_HVG, nrow(logged))
  top_genes <- order(gene_var, decreasing = TRUE)[seq_len(n_hvg)]
  hvg_mat <- as.matrix(logged[top_genes, , drop = FALSE])
  rm(logged); gc(verbose = FALSE)

  n_pcs <- min(N_PCS_REF, ncol(hvg_mat) - 1, nrow(hvg_mat) - 1)
  pca_fit <- irlba(t(hvg_mat), nv = n_pcs, center = TRUE, scale = FALSE)
  ref_embedding <- pca_fit$u %*% diag(pca_fit$d)
  rm(hvg_mat, pca_fit); gc(verbose = FALSE)

  true_knn <- BiocNeighbors::findKNN(ref_embedding, k = min(K_NEIGHBORS, nrow(ref_embedding) - 1))
  result <- list(true_knn_index = true_knn$index, n_cells = nrow(ref_embedding))
  saveRDS(result, cache_path)

  cat(sprintf("[%s] done in %s\n", inv$dataset_name[i], format(Sys.time() - t0)))
  rm(ref_embedding, true_knn, result); gc(verbose = FALSE)
}

cat("\nAll 6 real datasets' reference spaces cached.\n")
