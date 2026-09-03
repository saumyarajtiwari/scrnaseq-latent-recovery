# Step 6.5c v2 — Over-Smoothing detection, simulated graded-separability branch
#
# v1 (discarded) was mathematically circular: it projected group-mean
# transformed counts through each method's own fitted loadings, which is
# an unconditional matrix-algebra identity (group-mean and any fixed
# linear projection commute) regardless of whether the method preserves
# real signal -- explaining why pca_raw/libnorm/log all landed at ratio
# ~1.0000000 by construction, not by genuine measurement. It also used a
# wrong, assumed pca_shiftedlog transform formula (guessed log1p(x+1)
# rather than the actual data-driven per-file delta_opt shift), which
# would have needed fixing regardless, but the circularity made the whole
# approach unsalvageable independent of that error.
#
# v2 fixes this by comparing two genuinely independent quantities:
#   pre-PCA signal fraction:  ANOVA-style between-group SS / total SS,
#     computed directly on the method's transformed counts, in full gene
#     space (no group-averaging-then-projecting -- the source of the
#     circularity).
#   post-PCA signal fraction: identical formula, computed on the method's
#     actual embedding, using real per-cell true_group labels.
# Ratio = post-PCA fraction / pre-PCA fraction. No loadings, no
# re-derivation of any method's fitted transform is required for the
# post-PCA side (it's the method's real output); the two sides live in
# genuinely different spaces (raw gene expression vs. compressed PCs),
# so there is no shared-construction identity forcing the ratio toward 1.
#
# Scope: pca_raw, pca_libnorm, pca_log, pca_shiftedlog only (SCTransform v2
# and GLM-PCA still excluded -- their per-gene GLM-fitted residuals /
# nonlinear latent-factor models have no simple per-cell closed-form
# transform to apply to the full count matrix on the pre-PCA side).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
})

METHODS <- c("pca_raw", "pca_libnorm", "pca_log", "pca_shiftedlog")
TARGET_SUM <- 10000

log2fc <- read.csv("data/processed/step6_5b_simulated_log2fc.csv", stringsAsFactors = FALSE)
eligible <- log2fc[!is.na(log2fc$mean_abs_log2fc) &
                    (abs(log2fc$mean_abs_log2fc - 0.25) <= 0.05 |
                     abs(log2fc$mean_abs_log2fc - 0.50) <= 0.05), ]
eligible$target_band <- ifelse(abs(eligible$mean_abs_log2fc - 0.25) <= 0.05, "0.25", "0.50")
cat("Eligible unique raw files:", nrow(eligible), "\n")

# ANOVA-style signal fraction: sum over features (genes or PCs) of
# between-group SS, divided by sum over features of total SS. `mat` is
# features x cells (genes) or cells x features (embedding) -- handled via
# the `cells_are_rows` flag so the same function works for both spaces.
signal_fraction <- function(mat, groups, cells_are_rows) {
  if (!cells_are_rows) mat <- t(mat)  # normalize to cells x features
  grand_mean <- colMeans(mat)
  total_ss <- sum(sweep(mat, 2, grand_mean, "-")^2)
  if (total_ss == 0) return(NA_real_)
  between_ss <- 0
  for (g in unique(groups)) {
    idx <- groups == g
    n_g <- sum(idx)
    group_mean <- colMeans(mat[idx, , drop = FALSE])
    between_ss <- between_ss + n_g * sum((group_mean - grand_mean)^2)
  }
  between_ss / total_ss
}

derive_embedding_path <- function(raw_path, method) {
  bn <- basename(raw_path)
  sim <- basename(dirname(raw_path))
  suffix_map <- c(pca_raw = "rawpca", pca_libnorm = "libnormpca",
                   pca_log = "logpca", pca_shiftedlog = "shiftedlogpca")
  new_bn <- sub("\\.rds$", paste0("_", suffix_map[[method]], ".rds"), bn)
  file.path("data/processed", method, "simulated", sim, new_bn)
}

process_file <- function(i) {
  r <- eligible[i, ]
  raw_path <- r$file_path

  res <- tryCatch({
    src <- readRDS(raw_path)
    counts_mat <- counts(src)
    tg <- as.character(src$true_group)
    if (length(unique(tg)) < 2) stop("fewer than 2 groups")

    lib_sizes <- Matrix::colSums(counts_mat)
    safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
    scale_factor <- TARGET_SUM / safe_lib
    normalized <- counts_mat %*% Diagonal(x = scale_factor)

    delta_opt <- median(normalized@x)
    S <- normalized
    S@x <- log(S@x + delta_opt) - log(delta_opt)  # exact shiftedlog transform, per Step 3.4

    pre_pca_transforms <- list(
      pca_raw = counts_mat,
      pca_libnorm = normalized,
      pca_log = log1p(normalized),
      pca_shiftedlog = S
    )

    row_results <- list()
    for (method in METHODS) {
      m_res <- tryCatch({
        emb_path <- derive_embedding_path(raw_path, method)
        if (!file.exists(emb_path)) stop("embedding file not found: ", emb_path)
        emb_obj <- readRDS(emb_path)
        embedding <- emb_obj$embedding
        emb_tg <- as.character(emb_obj$true_group)

        tmat <- pre_pca_transforms[[method]]  # genes x cells
        pre_frac <- signal_fraction(as.matrix(tmat), tg, cells_are_rows = FALSE)
        post_frac <- signal_fraction(embedding, emb_tg, cells_are_rows = TRUE)

        ratio <- if (!is.na(pre_frac) && pre_frac > 0) post_frac / pre_frac else NA_real_

        list(pre_pca_signal_fraction = pre_frac, post_pca_signal_fraction = post_frac,
             ratio = ratio, over_smoothing_flag = isTRUE(ratio < 0.20),
             status = "ok", error_msg = NA_character_)
      }, error = function(e) {
        list(pre_pca_signal_fraction = NA, post_pca_signal_fraction = NA, ratio = NA,
             over_smoothing_flag = NA, status = "error", error_msg = conditionMessage(e))
      })
      m_res$method <- method
      row_results[[length(row_results) + 1]] <- m_res
    }
    row_results
  }, error = function(e) {
    list(list(method = NA, pre_pca_signal_fraction = NA, post_pca_signal_fraction = NA,
               ratio = NA, over_smoothing_flag = NA, status = "error",
               error_msg = paste("file-level error:", conditionMessage(e))))
  })

  lapply(res, function(x) { x$file_path <- raw_path; x$target_band <- r$target_band
                              x$mean_abs_log2fc <- r$mean_abs_log2fc; x })
}

cat("Starting mclapply (4 workers)...\n")
t0 <- Sys.time()
all_results <- mclapply(seq_len(nrow(eligible)), process_file, mc.cores = 4, mc.preschedule = FALSE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

flat_results <- do.call(rbind.data.frame, unlist(all_results, recursive = FALSE))
n_errors <- sum(flat_results$status == "error")
cat(sprintf("Total rows: %d | errors: %d\n", nrow(flat_results), n_errors))
if (n_errors > 0) print(table(flat_results$error_msg))

cat("\nRatio distribution by method:\n")
print(aggregate(ratio ~ method, data = flat_results, FUN = summary))
cat("\nOver-smoothing flagged (by method):\n")
print(aggregate(over_smoothing_flag ~ method, data = flat_results, FUN = mean))

write.csv(flat_results, "data/processed/step6_5c_simulated_variance_ratio.csv", row.names = FALSE)
cat("Written.\n")
