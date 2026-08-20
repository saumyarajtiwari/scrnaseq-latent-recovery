suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

RAW_PCA_NV <- 30   # was 50; fixed to match Step 3.9 spec (30 PCs, real data, cross-method comparability)
SEED       <- 42
OUT_DIR    <- "data/processed/pca_raw_real_fix"   # fresh location, never overwrite in place
LOG_FILE   <- file.path(OUT_DIR, "fix_progress.csv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_row <- function(file_path, status, msg = "") {
  row <- data.frame(timestamp = as.character(Sys.time()), file_path = file_path,
                     status = status, msg = msg, stringsAsFactors = FALSE)
  write.table(row, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

safe_saveRDS <- function(obj, out_file) {
  tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
  saveRDS(obj, tmp_file)
  file.rename(tmp_file, out_file)
}

real_files <- list(
  baron               = "data/real/baron/baron_harmonized.rds",
  muraro              = "data/real/muraro/muraro_harmonized.rds",
  pbmc68k             = "data/real/pbmc68k/pbmc68k_harmonized.rds",
  segerstolpe         = "data/real/segerstolpe/segerstolpe_harmonized.rds",
  tabula_sapiens_lung = "data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds",
  tasic2018           = "data/real/tasic2018/tasic2018_harmonized.rds"
)

fix_one <- function(name, source_file) {
  t0 <- Sys.time()
  out_dir <- file.path(OUT_DIR, name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(name, "_harmonized_rawpca.rds"))

  res <- tryCatch({
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
      embedding = embedding,
      singular_values = pca_fit$d,
      var_explained = var_explained,
      nv_used = nv,
      seed = SEED,
      method = "irlba_raw_pca",
      source_file = source_file,
      true_group = sce$true_group,
      batch_id = sce$batch_id,
      sim_params = NULL,
      loadings = loadings,
      loadings_method = "native_from_this_run",
      loadings_note = paste(
        "Regenerated 2026-08-20 to fix two issues in the original Step 3.1 real-data output:",
        "(1) nv was 50, inconsistent with Step 3.9's 30-PC real-data spec used by all other",
        "five methods; (2) loadings were never computed for real data in the original run or",
        "the subsequent Steps 3.1-3.4 backfill (which only covered simulated-grid files).",
        "This run recomputes embedding + loadings natively together via a single fresh irlba() call."
      ),
      flags = list(regenerated_reason = "nv_mismatch_and_missing_loadings_fix_2026-08-20"),
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    safe_saveRDS(out, out_file)
    list(status = "done", nv_used = nv, runtime = out$runtime_sec)
  }, error = function(e) list(status = "error", msg = conditionMessage(e)))

  if (res$status == "done") {
    log_row(out_file, "done", paste("nv_used:", res$nv_used, "runtime_sec:", round(res$runtime, 2)))
    cat("OK:", name, "| nv_used:", res$nv_used, "| runtime:", round(res$runtime, 2), "sec\n")
  } else {
    log_row(out_file, "error", res$msg)
    cat("ERROR:", name, "-", res$msg, "\n")
  }
}

cat("Regenerating 6 real-data raw-PCA files (nv=30, with native loadings)\n\n")
for (name in names(real_files)) {
  fix_one(name, real_files[[name]])
}
cat("\n=== DONE ===\n")
