suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

RAW_PCA_NV <- 50
SEED <- 42
LOG_FILE <- "data/processed/raw_recovery_progress.csv"

log_row <- function(file_path, status, msg = "") {
  row <- data.frame(timestamp = as.character(Sys.time()), file_path = file_path,
                     status = status, msg = msg, stringsAsFactors = FALSE)
  write.table(row, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

# Safe write: write to temp file in the SAME directory (ensures same filesystem
# for atomic rename), then rename over the original only on success
safe_saveRDS <- function(obj, out_file) {
  tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
  saveRDS(obj, tmp_file)
  file.rename(tmp_file, out_file)
}

recover_one <- function(out_file) {
  res <- tryCatch({
    # Reconstruct source file path from output filename convention
    fname <- basename(out_file)
    fname <- sub("_rawpca\\.rds$", "", fname)
    sim_name <- basename(dirname(out_file))
    source_file <- file.path("data/simulated/sce", sim_name, paste0(fname, ".rds"))
    if (!file.exists(source_file)) stop("source file not found: ", source_file)

    sce <- readRDS(source_file)
    counts_mat <- counts(sce)

    n_cells <- ncol(counts_mat); n_genes <- nrow(counts_mat)
    nv <- min(RAW_PCA_NV, n_cells - 1, n_genes - 1)

    set.seed(SEED)
    pca_fit <- irlba(t(counts_mat), nv = nv, center = FALSE, scale = FALSE)

    embedding <- pca_fit$u %*% diag(pca_fit$d)
    rownames(embedding) <- colnames(counts_mat)
    loadings <- pca_fit$v
    rownames(loadings) <- rownames(counts_mat)
    var_explained <- (pca_fit$d^2) / sum(pca_fit$d^2)

    out <- list(
      embedding = embedding, loadings = loadings, singular_values = pca_fit$d,
      var_explained = var_explained, nv_used = nv, seed = SEED,
      method = "irlba_raw_pca", source_file = source_file,
      true_group = sce$true_group, batch_id = sce$batch_id,
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      flags = list(note = "regenerated after corruption incident 2026-08-11; loadings native from this run, not backfilled"),
      recovery_timestamp = as.character(Sys.time())
    )
    safe_saveRDS(out, out_file)
    "done"
  }, error = function(e) paste("error:", conditionMessage(e)))

  if (res == "done") log_row(out_file, "done") else log_row(out_file, "error", res)
}

corrupted_files <- readRDS("data/corrupted_raw_files_list.rds")
cat("Regenerating", length(corrupted_files), "corrupted files\n")

for (i in seq_along(corrupted_files)) {
  recover_one(corrupted_files[i])
  if (i %% 500 == 0) cat("  progress:", i, "/", length(corrupted_files), "\n")
}

cat("\n=== RAW RECOVERY COMPLETE ===\n")
