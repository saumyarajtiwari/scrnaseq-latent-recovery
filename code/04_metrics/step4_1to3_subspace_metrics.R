suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(parallel)
})

METHOD_MAP <- list(
  pca_raw            = "rawpca",
  pca_libnorm        = "libnormpca",
  pca_log            = "logpca",
  pca_shiftedlog     = "shiftedlogpca",
  pca_sctransform_v2 = "sctpca",
  pca_glmpca         = "glmpca"
)

# Item 1 fix: ground-truth-linearity stratification, defined at the source
# where these scores are first computed (see PROJECT_HANDOVER.md ground-truth
# section). grassmann_distance/subspace_recovery_score/spectral_recovery_score
# are all computed against a LINEAR ground truth (tgm_centered's SVD basis)
# for every method. For the two nonlinear methods this partly measures "how
# linear is this method's recovered subspace," not purely biological
# recovery. No rows are removed; every downstream file that reads this
# script's output inherits method_family/ground_truth_caveat automatically.
METHOD_FAMILY <- c(
  pca_raw = "linear", pca_libnorm = "linear", pca_log = "linear",
  pca_shiftedlog = "linear", pca_sctransform_v2 = "nonlinear", pca_glmpca = "nonlinear"
)
GT_LINEARITY_CAVEAT <- paste(
  "subspace_recovery_score/grassmann_distance/spectral_recovery_score computed",
  "against a linearly-defined ground truth; for this nonlinear method the",
  "score partly reflects representational mismatch, not purely biological",
  "recovery -- see PROJECT_HANDOVER.md"
)

N_WORKERS  <- 4
BATCH_SIZE <- 500
OUT_DIR    <- "data/processed"
RESULTS_CSV    <- file.path(OUT_DIR, "step4_1to3_subspace_metrics.csv")
CHECKPOINT_RDS <- file.path(OUT_DIR, "step4_1to3_checkpoint.rds")
PROGRESS_LOG   <- file.path(OUT_DIR, "step4_1to3_progress.log")

log_msg <- function(...) {
  msg <- paste(format(Sys.time()), "-", ...)
  cat(msg, "\n")
  write(msg, PROGRESS_LOG, append = file.exists(PROGRESS_LOG))
}

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
sim_manifest <- manifest[manifest$data_type == "simulated" & !(manifest$is_null_control %in% TRUE), ]  # fix: == FALSE silently keeps NA rows as all-NA (same bug class Step 6.1 fixed); affects the 180 unparsed null-control-replicate rows (run_id=NA) discovered during Item 5's manifest rebuild

groups <- unique(sim_manifest[, c("source", "run_id")])
log_msg("Total (simulator, run_id) groups to process:", nrow(groups))

principal_angles_deg <- function(true_basis, est_basis) {
  Qt <- qr.Q(qr(true_basis))
  Qe <- qr.Q(qr(est_basis))
  s <- svd(t(Qt) %*% Qe)$d
  s <- pmin(pmax(s, -1), 1)
  acos(s) * 180 / pi
}

process_group <- function(sim_name, rid) {
  out_rows <- list()
  sce_path <- sprintf("data/simulated/sce/%s/%s_sce_run_%05d.rds", sim_name, sim_name, rid)
  res <- tryCatch({
    sce <- readRDS(sce_path)
    tgm_full <- metadata(sce)$true_group_means
    K <- ncol(tgm_full)
    r <- K - 1

    for (method in names(METHOD_MAP)) {
      methodtag <- METHOD_MAP[[method]]
      row_match <- sim_manifest[sim_manifest$source == sim_name & sim_manifest$run_id == rid & sim_manifest$method == method, ]
      if (nrow(row_match) != 1) {
        out_rows[[length(out_rows) + 1]] <- data.frame(
          source = sim_name, run_id = rid, method = method,
          method_family = METHOD_FAMILY[[method]],
          status = "error", error_msg = "manifest row not found",
          grassmann_distance = NA, subspace_recovery_score = NA,
          spectral_recovery_score = NA, max_principal_angle_deg = NA,
          n_flagged_slippage = NA,
          ground_truth_caveat = if (METHOD_FAMILY[[method]] == "nonlinear") GT_LINEARITY_CAVEAT else NA_character_,
          stringsAsFactors = FALSE
        )
        next
      }
      f_path <- row_match$file_path[1]

      m_res <- tryCatch({
        f <- readRDS(f_path)
        est_genes <- rownames(f$loadings)
        tgm_sub <- tgm_full[rownames(tgm_full) %in% est_genes, , drop = FALSE]
        tgm_sub <- tgm_sub[match(est_genes, rownames(tgm_sub)), , drop = FALSE]
        tgm_centered <- tgm_sub - rowMeans(tgm_sub)

        s_true <- svd(tgm_centered)
        true_basis <- s_true$u[, 1:r, drop = FALSE]
        true_dvals <- s_true$d[1:r]

        est_basis <- f$loadings[, 1:r, drop = FALSE]

        angles <- principal_angles_deg(true_basis, est_basis)
        angles_rad <- angles * pi / 180
        grassmann_dist <- sqrt(sum(sin(angles_rad)^2))
        subspace_score <- mean(cos(angles_rad)^2)
        spectral_score <- sum((true_dvals^2) * cos(angles_rad)^2) / sum(true_dvals^2)

        list(
          status = "ok", error_msg = NA,
          grassmann_distance = grassmann_dist,
          subspace_recovery_score = subspace_score,
          spectral_recovery_score = spectral_score,
          max_principal_angle_deg = max(angles),
          n_flagged_slippage = sum(angles > 30)
        )
      }, error = function(e) list(status = "error", error_msg = conditionMessage(e),
                                   grassmann_distance = NA, subspace_recovery_score = NA,
                                   spectral_recovery_score = NA, max_principal_angle_deg = NA,
                                   n_flagged_slippage = NA))

      out_rows[[length(out_rows) + 1]] <- data.frame(
        source = sim_name, run_id = rid, method = method,
        method_family = METHOD_FAMILY[[method]],
        status = m_res$status, error_msg = m_res$error_msg,
        grassmann_distance = m_res$grassmann_distance,
        subspace_recovery_score = m_res$subspace_recovery_score,
        spectral_recovery_score = m_res$spectral_recovery_score,
        max_principal_angle_deg = m_res$max_principal_angle_deg,
        n_flagged_slippage = m_res$n_flagged_slippage,
        ground_truth_caveat = if (METHOD_FAMILY[[method]] == "nonlinear") GT_LINEARITY_CAVEAT else NA_character_,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, out_rows)
  }, error = function(e) {
    data.frame(source = sim_name, run_id = rid, method = NA,
               method_family = NA_character_,
               status = "error", error_msg = paste("group-level:", conditionMessage(e)),
               grassmann_distance = NA, subspace_recovery_score = NA,
               spectral_recovery_score = NA, max_principal_angle_deg = NA,
               n_flagged_slippage = NA,
               ground_truth_caveat = NA_character_, stringsAsFactors = FALSE)
  })
  res
}

n_total <- nrow(groups)
n_batches <- ceiling(n_total / BATCH_SIZE)
completed_batches <- if (file.exists(CHECKPOINT_RDS)) readRDS(CHECKPOINT_RDS) else integer(0)
log_msg("Resuming:", length(completed_batches), "of", n_batches, "batches already done")

for (b in seq_len(n_batches)) {
  if (b %in% completed_batches) next
  idx_start <- (b - 1) * BATCH_SIZE + 1
  idx_end <- min(b * BATCH_SIZE, n_total)
  batch_groups <- groups[idx_start:idx_end, ]

  results <- mclapply(seq_len(nrow(batch_groups)), function(i) {
    process_group(batch_groups$source[i], batch_groups$run_id[i])
  }, mc.cores = N_WORKERS, mc.preschedule = FALSE)

  batch_df <- do.call(rbind, results)
  write.table(batch_df, RESULTS_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(RESULTS_CSV), append = file.exists(RESULTS_CSV))

  completed_batches <- c(completed_batches, b)
  saveRDS(completed_batches, CHECKPOINT_RDS)

  n_errors <- sum(batch_df$status == "error", na.rm = TRUE)
  log_msg(sprintf("Batch %d/%d done (%d groups, %d rows). Errors: %d. Total groups processed: %d/%d",
                   b, n_batches, nrow(batch_groups), nrow(batch_df), n_errors, idx_end, n_total))
}

log_msg("=== STEP 4.1-4.3 COMPLETE ===")
final <- read.csv(RESULTS_CSV, stringsAsFactors = FALSE)
log_msg("Total rows:", nrow(final), "| Total errors:", sum(final$status == "error"))
