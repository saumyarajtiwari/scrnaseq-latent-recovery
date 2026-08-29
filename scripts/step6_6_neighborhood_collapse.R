# Step 6.6 — Detect Neighborhood Collapse (v3: atomic cache writes)
#
# [... same header as before ...]
#
# BUG FIX 2 (found after full run: 963 errors, 0.49%, clustering by shared
# run_id across multiple methods -- e.g. run_id 3143 (splatter): pca_raw,
# pca_sctransform_v2, pca_glmpca, pca_libnorm all failed together; run_id
# 8871 (scdesign3): all 4 remaining methods failed together. This is the
# signature of a shared, torn cache file: build_reference_space() used a
# plain saveRDS() with no atomic temp-file-then-rename, violating this
# project's established "atomic writes are mandatory for any script
# modifying/creating files" convention (from the Step 1/3 disk-full
# corruption incidents). Multiple workers processing different methods for
# the SAME source file could each see file.exists()==FALSE simultaneously
# (TOCTOU race) and write to the identical cache path concurrently, leaving
# a torn file that fails for every subsequent reader -- consistent with the
# "ReadItem: unknown type X" and "error reading from connection" errors
# observed, and with errors clustering by shared run_id rather than being
# randomly distributed. Fixed via temp-file + atomic rename: concurrent
# workers may still redundantly rebuild the same reference space, but the
# final file is now always a complete, valid write from exactly one worker.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(BiocNeighbors)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args
retry_only <- "--retry_errors" %in% args
set.seed(42)

K_NEIGHBORS <- 15
N_HVG <- 500
N_PCS_REF <- 15
OVERLAP_THRESHOLD <- 0.5

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
work <- manifest[!(manifest$is_null_control %in% TRUE) & manifest$status == "ok", ]

master <- read.csv("data/processed/step4_master_results_table.csv", stringsAsFactors = FALSE)
master_key <- paste(master$data_type, master$source, master$run_id, master$method)
cont_lookup <- setNames(master$continuity, master_key)

cat(sprintf("Total candidate rows: %d\n", nrow(work)))

if (retry_only) {
  prior <- read.csv("data/processed/step6_6_neighborhood_collapse.csv", stringsAsFactors = FALSE)
  err_paths <- prior$file_path[prior$status == "error"]
  work <- work[work$file_path %in% err_paths, ]
  cat(sprintf("RETRY MODE: %d previously-errored rows\n", nrow(work)))
} else if (dry_run) {
  set.seed(1)
  n_sample <- 500
  large_real <- work[work$source %in% c("pbmc68k", "tabula_sapiens_lung", "tasic2018"), ]
  n_large <- min(9, nrow(large_real))
  large_real_sample <- large_real[sample(nrow(large_real), n_large), ]
  remaining <- work[!(rownames(work) %in% rownames(large_real_sample)), ]
  rest_sample <- remaining[sample(nrow(remaining), n_sample - n_large), ]
  work <- rbind(large_real_sample, rest_sample)
  cat(sprintf("DRY RUN: subset to %d rows (incl. %d large-real rows)\n", nrow(work), n_large))
}

ref_cache_dir <- "data/processed/step6_reference_space_cache"
dir.create(ref_cache_dir, showWarnings = FALSE, recursive = TRUE)

build_reference_space <- function(source_file) {
  cache_key <- gsub("[/\\\\]", "_", source_file)
  cache_path <- file.path(ref_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(cache_path)) {
    result <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(result)) return(result)
    # Cache file exists but is corrupted (torn write from prior race
    # condition) -- fall through and rebuild it properly below.
  }

  src <- readRDS(source_file)
  counts_mat <- counts(src)
  lib_sizes <- Matrix::colSums(counts_mat)
  safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
  normalized <- counts_mat %*% Diagonal(x = 10000 / safe_lib)
  logged <- log1p(normalized)

  n_cells_local <- ncol(logged)
  gene_mean <- Matrix::rowMeans(logged)
  gene_mean_sq <- Matrix::rowMeans(logged^2)
  gene_var <- (gene_mean_sq - gene_mean^2) * (n_cells_local / (n_cells_local - 1))
  n_hvg <- min(N_HVG, nrow(logged))
  top_genes <- order(gene_var, decreasing = TRUE)[seq_len(n_hvg)]
  hvg_mat <- as.matrix(logged[top_genes, , drop = FALSE])

  n_pcs <- min(N_PCS_REF, ncol(hvg_mat) - 1, nrow(hvg_mat) - 1)
  pca_fit <- irlba::irlba(t(hvg_mat), nv = n_pcs, center = TRUE, scale = FALSE)
  ref_embedding <- pca_fit$u %*% diag(pca_fit$d)

  true_knn <- BiocNeighbors::findKNN(ref_embedding, k = min(K_NEIGHBORS, nrow(ref_embedding) - 1))

  result <- list(true_knn_index = true_knn$index, n_cells = nrow(ref_embedding))

  # ATOMIC WRITE: temp file (unique per-process via PID) + rename, so
  # concurrent workers can never produce a torn/partial cache file.
  tmp_path <- paste0(cache_path, ".tmp_", Sys.getpid())
  saveRDS(result, tmp_path)
  file.rename(tmp_path, cache_path)
  result
}

process_row <- function(i) {
  r <- work[i, ]
  key <- paste(r$data_type, r$source, r$run_id, r$method)
  cont <- cont_lookup[[key]]

  out <- list(
    data_type = r$data_type, source = r$source, run_id = r$run_id, method = r$method,
    file_path = r$file_path, n_cells = NA, continuity = ifelse(is.null(cont), NA, cont),
    continuity_trigger = NA, mean_overlap_score = NA, overlap_flag = NA,
    status = "ok", error_msg = NA_character_
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    n_total <- nrow(emb)

    ref <- build_reference_space(obj$source_file)
    if (ref$n_cells != n_total) stop(sprintf("alignment failure: reference has %d cells, embedding has %d", ref$n_cells, n_total))

    k_use <- min(K_NEIGHBORS, n_total - 1)
    emb_knn <- BiocNeighbors::findKNN(emb, k = k_use)

    overlaps <- sapply(seq_len(n_total), function(cell) {
      true_set <- ref$true_knn_index[cell, ]
      emb_set <- emb_knn$index[cell, ]
      length(intersect(true_set, emb_set)) / k_use
    })
    mean_overlap <- mean(overlaps)

    list(n_cells = n_total, continuity_trigger = isTRUE(cont < 0.85),
         mean_overlap_score = mean_overlap, overlap_flag = isTRUE(mean_overlap < OVERLAP_THRESHOLD),
         status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(n_cells = NA, continuity_trigger = isTRUE(cont < 0.85), mean_overlap_score = NA,
         overlap_flag = NA, status = "error", error_msg = conditionMessage(e))
  })

  out[names(res)] <- res
  out
}

cat("Starting mclapply (4 workers, mc.preschedule=FALSE)...\n")
t0 <- Sys.time()
results_list <- mclapply(seq_len(nrow(work)), process_row, mc.cores = 4, mc.preschedule = FALSE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

results_df <- do.call(rbind.data.frame, results_list)
n_errors <- sum(results_df$status == "error")
cat(sprintf("Rows processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) print(table(results_df$error_msg))

if (retry_only) {
  prior <- read.csv("data/processed/step6_6_neighborhood_collapse.csv", stringsAsFactors = FALSE)
  key_cols <- c("data_type","source","run_id","method","file_path")
  prior_key <- do.call(paste, prior[key_cols])
  retry_key <- do.call(paste, results_df[key_cols])
  match_idx <- match(retry_key, prior_key)
  prior[match_idx, names(results_df)] <- results_df
  write.csv(prior, "data/processed/step6_6_neighborhood_collapse.csv", row.names = FALSE)
  cat("Merged retry results into main output. Final error count:",
      sum(prior$status == "error"), "\n")
} else {
  cat("\nContinuity-triggered:", sum(results_df$continuity_trigger, na.rm=TRUE), "/", nrow(results_df), "\n")
  cat("Overlap-flagged:", sum(results_df$overlap_flag, na.rm=TRUE), "/", nrow(results_df), "\n")
  suffix <- if (dry_run) "_DRYRUN" else ""
  write.csv(results_df, paste0("data/processed/step6_6_neighborhood_collapse", suffix, ".csv"), row.names = FALSE)
  cat("Written.\n")
}
