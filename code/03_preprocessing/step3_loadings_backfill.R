suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

TARGET_SUM <- 10000
LOG_FILE <- "data/processed/loadings_backfill_progress.csv"
DISK_MOUNT <- "/mnt/extra2"
MIN_FREE_GB <- 2.0
CHECK_DISK_EVERY <- 1000

log_row <- function(method, file_path, status, msg = "") {
  row <- data.frame(timestamp = as.character(Sys.time()), method = method,
                     file_path = file_path, status = status, msg = msg,
                     stringsAsFactors = FALSE)
  write.table(row, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

get_free_gb <- function() {
  out <- system(paste("df --output=avail", DISK_MOUNT, "| tail -1"), intern = TRUE)
  as.numeric(trimws(out)) / 1e6
}

# Safe write: temp file in same dir, atomic rename on success
safe_saveRDS <- function(obj, out_file) {
  tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
  saveRDS(obj, tmp_file)
  file.rename(tmp_file, out_file)
}

get_transformed_matrix <- function(method, sce) {
  counts_mat <- counts(sce)
  lib_sizes <- Matrix::colSums(counts_mat)
  safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
  scale_factor <- TARGET_SUM / safe_lib
  normalized <- counts_mat %*% Diagonal(x = scale_factor)
  if (method == "libnorm") return(t(normalized))
  else if (method == "log") return(t(log1p(normalized)))
}

backfill_one <- function(method, out_file, centered) {
  obj <- tryCatch(readRDS(out_file), error = function(e) NULL)
  if (is.null(obj)) { gc(verbose = FALSE); log_row(method, out_file, "error", "could not read output file"); return(invisible()) }
  if (!is.null(obj$loadings)) { rm(obj); gc(verbose = FALSE); log_row(method, out_file, "skipped_existing"); return(invisible()) }

  res <- tryCatch({
    sce <- readRDS(obj$source_file)

    if (method == "shiftedlog") {
      counts_mat <- counts(sce)
      lib_sizes <- Matrix::colSums(counts_mat)
      safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
      normalized <- counts_mat %*% Diagonal(x = TARGET_SUM / safe_lib)
      delta_opt <- obj$delta_used
      S <- normalized
      S@x <- log(S@x + delta_opt) - log(delta_opt)
      A <- t(S)
    } else {
      A <- get_transformed_matrix(method, sce)
    }
    d <- obj$singular_values
    embedding <- obj$embedding
    col_means <- Matrix::colMeans(A)
    if (centered) {
      V <- (Matrix::crossprod(A, embedding) - outer(col_means, Matrix::colSums(embedding))) %*% diag(1 / d^2)
    } else {
      V <- Matrix::crossprod(A, embedding) %*% diag(1 / d^2)
    }
    rownames(V) <- rownames(counts(sce))

    obj$loadings <- V
    obj$loadings_method <- "analytically_reconstructed_from_embedding"
    obj$loadings_note <- paste(
      "Loadings derived post-hoc via V = A_centered^T x embedding x diag(1/d^2),",
      "exact under SVD identities. Precision decreases for smaller singular values;",
      "see stored singular_values to judge per-component reliability."
    )
    safe_saveRDS(obj, out_file)
    rm(sce, A, embedding, V, obj); gc(verbose = FALSE)
    "done"
  }, error = function(e) paste("error:", conditionMessage(e)))

  if (res == "done") log_row(method, out_file, "done") else log_row(method, out_file, "error", res)
}

methods <- list(
  libnorm    = list(dir_sim = "data/processed/pca_libnorm/simulated", dir_real = "data/processed/pca_libnorm/real", centered = TRUE),
  log        = list(dir_sim = "data/processed/pca_log/simulated", dir_real = "data/processed/pca_log/real", centered = TRUE),
  shiftedlog = list(dir_sim = "data/processed/pca_shiftedlog/simulated", dir_real = "data/processed/pca_shiftedlog/real", centered = TRUE)
)

for (m in names(methods)) {
  cfg <- methods[[m]]
  files <- c(
    list.files(cfg$dir_sim, recursive = TRUE, full.names = TRUE, pattern = "\\.rds$"),
    list.files(cfg$dir_real, recursive = TRUE, full.names = TRUE, pattern = "\\.rds$")
  )
  cat("=== Backfilling", m, "(", length(files), "files ) ===\n")
  for (i in seq_along(files)) {
    backfill_one(m, files[i], cfg$centered)
    if (i %% 500 == 0) gc(verbose = FALSE)
    if (i %% CHECK_DISK_EVERY == 0) {
      fg <- get_free_gb()
      cat("  ", m, "progress:", i, "/", length(files), " free disk (GB):", round(fg, 2), "\n")
      if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free on ", DISK_MOUNT, ". Aborting.")
    }
  }
  cat(m, "complete.\n\n")
}

cat("=== LOADINGS BACKFILL COMPLETE ===\n")
