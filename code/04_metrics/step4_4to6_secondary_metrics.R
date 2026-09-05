suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(FNN)
  library(mclust)
  library(irlba)
  library(parallel)
})

METHOD_MAP <- list(
  pca_raw = "rawpca", pca_libnorm = "libnormpca", pca_log = "logpca",
  pca_shiftedlog = "shiftedlogpca", pca_sctransform_v2 = "sctpca", pca_glmpca = "glmpca"
)

K_NEIGHBORS   <- 15
K_EXT         <- 100
N_COMPONENTS  <- 30
REF_N_HVG     <- 500
REF_N_PCS     <- 15
SEED          <- 42
N_WORKERS_SIM  <- 12
N_WORKERS_REAL <- 1
BATCH_SIZE_SIM <- 1000
OUT_DIR <- "data/processed"
RESULTS_CSV    <- file.path(OUT_DIR, "step4_4to6_secondary_metrics.csv")
CHECKPOINT_SIM <- file.path(OUT_DIR, "step4_4to6_checkpoint_sim.rds")
CHECKPOINT_REAL <- file.path(OUT_DIR, "step4_4to6_checkpoint_real.rds")
PROGRESS_LOG   <- file.path(OUT_DIR, "step4_4to6_progress.log")

log_msg <- function(...) {
  msg <- paste(format(Sys.time()), "-", ...)
  cat(msg, "\n")
  write(msg, PROGRESS_LOG, append = file.exists(PROGRESS_LOG))
}

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)

compute_reference <- function(counts_mat, n_hvg = REF_N_HVG, n_pcs = REF_N_PCS, n_prefilter = 3000) {
  # Two-stage gene reduction for datasets with unusually high gene counts
  # (Tabula Sapiens Lung: 56,139; Tasic2018: 42,865, vs PBMC68k's 20,387).
  # A cheap mean-expression pre-selection (top 3000, same cap already
  # established and validated for SCTransform/GLM-PCA in Step 3.5/3.6 on
  # these same two datasets) happens BEFORE any log-transform/variance
  # computation, since holding multiple full-size sparse copies at the
  # original gene count was confirmed via dmesg to OOM-kill even alone,
  # single-threaded, and a lenient >=3-cell detection filter alone was
  # insufficient (only removed ~10% of genes for these dense datasets).
  if (nrow(counts_mat) > n_prefilter) {
    gene_mean_counts <- Matrix::rowMeans(counts_mat)
    top_idx <- order(gene_mean_counts, decreasing = TRUE)[1:n_prefilter]
    counts_mat <- counts_mat[top_idx, , drop = FALSE]
    gc(verbose = FALSE)
  }

  lib_sizes <- Matrix::colSums(counts_mat); lib_sizes[lib_sizes == 0] <- 1
  norm_mat <- Matrix::t(Matrix::t(counts_mat) / lib_sizes) * 10000
  log_mat <- log1p(norm_mat)
  rm(norm_mat); gc(verbose = FALSE)
  gene_means <- Matrix::rowMeans(log_mat)
  gene_vars <- Matrix::rowMeans(log_mat^2) - gene_means^2
  n_hvg_eff <- min(n_hvg, nrow(log_mat))
  hvg_idx <- order(gene_vars, decreasing = TRUE)[1:n_hvg_eff]
  log_mat_hvg <- log_mat[hvg_idx, , drop = FALSE]
  rm(log_mat); gc(verbose = FALSE)
  n_pcs_eff <- min(n_pcs, ncol(log_mat_hvg) - 1, nrow(log_mat_hvg) - 1)
  set.seed(SEED)
  pca_res <- irlba(t(log_mat_hvg), nv = n_pcs_eff, center = TRUE, scale = FALSE)
  pca_res$u %*% diag(pca_res$d)
}

trust_cont <- function(nnX_idx, nnY_idx, k, k_ext, n) {
  T_sum <- 0; C_sum <- 0
  for (i in 1:n) {
    NkX <- nnX_idx[i, 1:k]; NkY <- nnY_idx[i, 1:k]
    intruders <- setdiff(NkY, NkX)
    for (j in intruders) {
      pos <- match(j, nnX_idx[i, ]); rk <- if (is.na(pos)) k_ext + 1 else pos
      T_sum <- T_sum + (rk - k)
    }
    missing_pts <- setdiff(NkX, NkY)
    for (j in missing_pts) {
      pos <- match(j, nnY_idx[i, ]); rk <- if (is.na(pos)) k_ext + 1 else pos
      C_sum <- C_sum + (rk - k)
    }
  }
  norm_const <- 2 / (n * k * (2 * n - 3 * k - 1))
  list(trustworthiness = 1 - norm_const * T_sum, continuity = 1 - norm_const * C_sum)
}

simplified_silhouette <- function(embedding, labels) {
  groups <- unique(labels)
  if (length(groups) < 2) return(NA_real_)
  centroids <- t(sapply(groups, function(g) colMeans(embedding[labels == g, , drop = FALSE])))
  x2 <- rowSums(embedding^2); c2 <- rowSums(centroids^2)
  dist_mat <- sqrt(pmax(outer(x2, c2, "+") - 2 * embedding %*% t(centroids), 0))
  own_idx <- match(labels, groups)
  a_i <- dist_mat[cbind(seq_along(labels), own_idx)]
  dist_mat_others <- dist_mat
  dist_mat_others[cbind(seq_along(labels), own_idx)] <- Inf
  b_i <- apply(dist_mat_others, 1, min)
  mean((b_i - a_i) / pmax(a_i, b_i), na.rm = TRUE)
}

process_group <- function(data_type, source_name, run_id, sce_path, method_rows) {
  out_rows <- list()
  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)
    true_group_full <- sce$true_group
    n <- ncol(counts_mat)
    ref <- compute_reference(counts_mat)
    k_ext_eff <- min(K_EXT, n - 1)
    k_use <- min(K_NEIGHBORS, k_ext_eff)
    set.seed(SEED)
    nnX <- get.knn(ref, k = k_ext_eff)$nn.index

    has_label <- !is.na(true_group_full)
    K_true <- length(unique(true_group_full[has_label]))

    for (method in names(METHOD_MAP)) {
      f_path <- method_rows$file_path[method_rows$method == method]
      m_res <- tryCatch({
        if (length(f_path) != 1) stop("manifest row not found")
        f <- readRDS(f_path)
        emb <- f$embedding[, 1:min(N_COMPONENTS, ncol(f$embedding)), drop = FALSE]

        nnY <- get.knn(emb, k = k_ext_eff)$nn.index
        tc <- trust_cont(nnX, nnY, k_use, k_ext_eff, n)

        ari <- NA_real_; sil <- NA_real_; label_note <- "ok"
        if (K_true >= 2) {
          emb_labeled <- emb[has_label, , drop = FALSE]
          labels_labeled <- true_group_full[has_label]
          set.seed(SEED)
          km <- kmeans(emb_labeled, centers = K_true, nstart = 10)
          ari <- adjustedRandIndex(km$cluster, labels_labeled)
          sil <- simplified_silhouette(emb_labeled, labels_labeled)
        } else {
          label_note <- "no_labels"
        }

        list(status = "ok", error_msg = NA, label_note = label_note,
             trustworthiness = tc$trustworthiness, continuity = tc$continuity,
             ari = ari, silhouette = sil, n_labeled = sum(has_label))
      }, error = function(e) list(status = "error", error_msg = conditionMessage(e), label_note = NA,
                                   trustworthiness = NA, continuity = NA, ari = NA, silhouette = NA, n_labeled = NA))

      out_rows[[length(out_rows) + 1]] <- data.frame(
        data_type = data_type, source = source_name, run_id = run_id, method = method,
        status = m_res$status, error_msg = m_res$error_msg, label_note = m_res$label_note,
        trustworthiness = m_res$trustworthiness, continuity = m_res$continuity,
        ari = m_res$ari, silhouette = m_res$silhouette, n_labeled = m_res$n_labeled,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, out_rows)
  }, error = function(e) {
    data.frame(data_type = data_type, source = source_name, run_id = run_id, method = NA,
               status = "error", error_msg = paste("group-level:", conditionMessage(e)), label_note = NA,
               trustworthiness = NA, continuity = NA, ari = NA, silhouette = NA, n_labeled = NA,
               stringsAsFactors = FALSE)
  })
  res
}

write_batch <- function(batch_df) {
  write.table(batch_df, RESULTS_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(RESULTS_CSV), append = file.exists(RESULTS_CSV))
}

sim_manifest <- manifest[manifest$data_type == "simulated" & !(manifest$is_null_control %in% TRUE), ]  # fix: == FALSE silently keeps NA rows as all-NA (same bug class Step 6.1 fixed); affects the 180 unparsed null-control-replicate rows (run_id=NA) discovered during Item 5's manifest rebuild
sim_groups <- unique(sim_manifest[, c("source", "run_id")])
sim_groups$sce_path <- sprintf("data/simulated/sce/%s/%s_sce_run_%05d.rds",
                                sim_groups$source, sim_groups$source, sim_groups$run_id)
log_msg("PHASE 1 (simulated): total groups:", nrow(sim_groups))

n_total_sim <- nrow(sim_groups)
n_batches_sim <- ceiling(n_total_sim / BATCH_SIZE_SIM)
completed_sim <- if (file.exists(CHECKPOINT_SIM)) readRDS(CHECKPOINT_SIM) else integer(0)
log_msg("Resuming simulated phase:", length(completed_sim), "of", n_batches_sim, "batches already done")

for (b in seq_len(n_batches_sim)) {
  if (b %in% completed_sim) next
  idx_start <- (b - 1) * BATCH_SIZE_SIM + 1
  idx_end <- min(b * BATCH_SIZE_SIM, n_total_sim)
  batch_groups <- sim_groups[idx_start:idx_end, ]

  results <- mclapply(seq_len(nrow(batch_groups)), function(i) {
    g <- batch_groups[i, ]
    method_rows <- sim_manifest[sim_manifest$source == g$source & sim_manifest$run_id == g$run_id, c("method","file_path")]
    process_group("simulated", g$source, g$run_id, g$sce_path, method_rows)
  }, mc.cores = N_WORKERS_SIM, mc.preschedule = FALSE)

  batch_df <- do.call(rbind, results)
  write_batch(batch_df)
  completed_sim <- c(completed_sim, b)
  saveRDS(completed_sim, CHECKPOINT_SIM)

  n_errors <- sum(batch_df$status == "error", na.rm = TRUE)
  log_msg(sprintf("[SIM] Batch %d/%d done (%d groups, %d rows). Errors: %d. Processed: %d/%d",
                   b, n_batches_sim, nrow(batch_groups), nrow(batch_df), n_errors, idx_end, n_total_sim))
}
log_msg("=== PHASE 1 (simulated) COMPLETE ===")

real_manifest <- manifest[manifest$data_type == "real", ]
real_sources <- unique(real_manifest$source)
log_msg("PHASE 2 (real): total datasets:", length(real_sources))

completed_real <- if (file.exists(CHECKPOINT_REAL)) readRDS(CHECKPOINT_REAL) else character(0)
remaining_real <- setdiff(real_sources, completed_real)
log_msg("Resuming real phase:", length(completed_real), "of", length(real_sources), "datasets already done")

if (length(remaining_real) > 0) {
  results <- mclapply(remaining_real, function(src) {
    method_rows <- real_manifest[real_manifest$source == src, c("method","file_path")]
    sce_path <- sprintf("data/real/%s/%s_harmonized.rds", src, src)
    process_group("real", src, NA_integer_, sce_path, method_rows)
  }, mc.cores = N_WORKERS_REAL, mc.preschedule = FALSE)

  for (i in seq_along(remaining_real)) {
    r <- results[[i]]
    is_valid <- is.data.frame(r) && nrow(r) == length(METHOD_MAP) && all(c("status","source") %in% names(r))
    if (is_valid == FALSE) {
      log_msg(sprintf("[REAL] %s FAILED - worker process likely killed (invalid/missing result, class=%s). NOT marked complete, will retry on next run.", remaining_real[i], paste(class(r), collapse=",")))
      next
    }
    write_batch(r)
    completed_real <- c(completed_real, remaining_real[i])
    saveRDS(completed_real, CHECKPOINT_REAL)
    n_errors <- sum(r$status == "error", na.rm = TRUE)
    log_msg(sprintf("[REAL] %s done (%d rows). Errors: %d.", remaining_real[i], nrow(r), n_errors))
  }
}
log_msg("=== PHASE 2 (real) COMPLETE ===")

log_msg("=== STEP 4.4-4.6 FULLY COMPLETE ===")
final <- read.csv(RESULTS_CSV, stringsAsFactors = FALSE)
log_msg("Total rows:", nrow(final), "| Total errors:", sum(final$status == "error"))
