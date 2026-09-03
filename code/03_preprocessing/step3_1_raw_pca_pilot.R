suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

RAW_PCA_NV <- 50
SEED <- 42

run_raw_pca <- function(sce_path, out_dir, flags = list()) {
  cat("=== ", sce_path, " ===\n")
  t0 <- Sys.time()

  sce <- readRDS(sce_path)
  counts_mat <- counts(sce)  # dgCMatrix, genes x cells

  s <- summary(counts_mat)
  is_integer_like <- all(s$x == round(s$x))
  cat("  dim (genes x cells):", paste(dim(counts_mat), collapse=" x "), "\n")
  cat("  nnzero:", length(s$x), " integer-valued (spot check):", is_integer_like, "\n")

  n_cells <- ncol(counts_mat)
  n_genes <- nrow(counts_mat)
  nv <- min(RAW_PCA_NV, n_cells - 1, n_genes - 1)

  set.seed(SEED)
  pca_fit <- irlba(t(counts_mat), nv = nv, center = FALSE, scale = FALSE)

  embedding <- pca_fit$u %*% diag(pca_fit$d)
  rownames(embedding) <- colnames(counts_mat)
  var_explained <- (pca_fit$d^2) / sum(pca_fit$d^2)

  out <- list(
    embedding = embedding,
    singular_values = pca_fit$d,
    var_explained = var_explained,
    nv_used = nv,
    seed = SEED,
    method = "irlba_raw_pca",
    source_file = sce_path,
    true_group = sce$true_group,
    batch_id = sce$batch_id,
    flags = flags,
    runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_rawpca.rds"))
  saveRDS(out, out_file)

  cat("  nv used:", nv, " runtime (s):", round(out$runtime_sec, 2),
      " output size (MB):", round(file.size(out_file)/1e6, 3), "\n")

  rm(sce, counts_mat, pca_fit, embedding)
  gc(verbose = FALSE)
  cat("\n")

  out_file
}

pilot_targets <- list(
  list(path = "data/simulated/sce/scdesign3/scdesign3_sce_run_00001.rds",
       out  = "data/processed/pca_raw/simulated/scdesign3"),
  list(path = list.files("data/simulated/sce/splatter", full.names = TRUE, pattern = "\\.rds$")[1],
       out  = "data/processed/pca_raw/simulated/splatter"),
  list(path = list.files("data/simulated/sce/symsim", full.names = TRUE, pattern = "\\.rds$")[1],
       out  = "data/processed/pca_raw/simulated/symsim"),
  list(path = "data/real/pbmc68k/pbmc68k_harmonized.rds",
       out  = "data/processed/pca_raw/real/pbmc68k")
)

results <- list()
for (target in pilot_targets) {
  results[[target$path]] <- tryCatch(
    run_raw_pca(target$path, target$out),
    error = function(e) {
      cat("  ERROR on", target$path, ":", conditionMessage(e), "\n\n")
      NA
    }
  )
}

cat("=== PILOT SUMMARY ===\n")
print(results)
