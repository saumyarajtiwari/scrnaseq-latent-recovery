# Step 6.3 v2 — Detect Phantom Clustering (with empirically-calibrated CH threshold)
#
# Revision rationale: the literal spec threshold (CH > 10) was tested against
# a pure-noise baseline (20 draws, n=1000, 30 dims) and found to sit well
# BELOW what pure i.i.d. Gaussian noise alone produces under forced k=3
# k-means (baseline CH ~26.6-27.8) -- CH has no built-in chance correction,
# unlike ARI, so an uncalibrated absolute threshold does not discriminate
# real artifactual structure from pure noise at this n/k. This is standard
# practice for un-normalized clustering-validity indices: calibrate against
# a matched-null (permutation/simulation) distribution rather than an
# absolute cutoff (e.g. Tibshirani et al.'s Gap Statistic, and standard
# permutation-test practice generally).
#
# Both flags are reported: phantom_flag_ch_literal (per spec, CH > 10) and
# phantom_flag_ch_relative (CH > empirical 95th percentile of matched-
# dimensionality Gaussian-noise baseline), the latter being the
# scientifically meaningful indicator.
#
# Also fixes: kmeans non-convergence warning seen in v1, by raising
# iter.max (default 10 was insufficient at this n/k in some cases).

suppressPackageStartupMessages(library(aricode))
set.seed(42)

nc_manifest <- read.csv("data/processed/step6_3_null_control_manifest.csv", stringsAsFactors = FALSE)
cat("Total null-control files:", nrow(nc_manifest), "\n")

calinski_harabasz <- function(emb, cluster) {
  k <- length(unique(cluster)); n <- nrow(emb)
  overall <- colMeans(emb)
  between_ss <- 0; within_ss <- 0
  for (c in unique(cluster)) {
    pts <- emb[cluster == c, , drop = FALSE]
    cen <- colMeans(pts)
    between_ss <- between_ss + nrow(pts) * sum((cen - overall)^2)
    within_ss <- within_ss + sum(sweep(pts, 2, cen, "-")^2)
  }
  if (within_ss == 0) return(Inf)
  (between_ss / (k - 1)) / (within_ss / (n - k))
}

## ---- Step 1: empirical noise-baseline calibration, per embedding dimensionality ----
cat("\n=== Calibrating noise baselines ===\n")
dims_needed <- c(30, 50)  # confirmed: pca_raw=50, all other 5 methods=30
N_BASELINE_DRAWS <- 200
noise_baseline <- list()

for (d in dims_needed) {
  set.seed(1000 + d)
  ch_vals <- replicate(N_BASELINE_DRAWS, {
    noise <- matrix(rnorm(1000 * d), nrow = 1000)
    km <- kmeans(noise, centers = 3, nstart = 25, iter.max = 100)$cluster
    calinski_harabasz(noise, km)
  })
  noise_baseline[[as.character(d)]] <- ch_vals
  cat(sprintf("dim=%d: noise CH mean=%.2f, 95th pct=%.2f, 99th pct=%.2f, max=%.2f\n",
              d, mean(ch_vals), quantile(ch_vals, 0.95), quantile(ch_vals, 0.99), max(ch_vals)))
}

## ---- Step 2: main computation on real null-control files ----
process_row <- function(i) {
  r <- nc_manifest[i, ]
  out <- list(
    method = r$method, source = r$source, run_id = r$run_id, replicate = r$replicate,
    sparsity = r$sparsity, n_cells = NA, embedding_dim = NA,
    ari_vs_random = NA, ch_index = NA, ch_noise_95pct = NA,
    phantom_flag_ari = NA, phantom_flag_ch_literal = NA, phantom_flag_ch_relative = NA,
    phantom_flag_any_relative = NA, status = "ok", error_msg = NA_character_
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    n <- nrow(emb); d <- ncol(emb)

    set.seed(42)
    km <- kmeans(emb, centers = 3, nstart = 25, iter.max = 100)$cluster

    set.seed(42 + r$run_id * 10 + r$replicate)
    random_labels <- sample(1:3, n, replace = TRUE)
    ari <- aricode::ARI(as.integer(km), as.integer(random_labels))

    ch <- calinski_harabasz(emb, km)
    baseline_key <- as.character(d)
    ch_95pct <- if (baseline_key %in% names(noise_baseline)) {
      quantile(noise_baseline[[baseline_key]], 0.95)
    } else {
      warning(sprintf("No noise baseline for dim=%d, file=%s", d, r$file_path))
      NA
    }

    flag_ari <- isTRUE(ari > 0.2)
    flag_ch_lit <- isTRUE(ch > 10)
    flag_ch_rel <- isTRUE(!is.na(ch_95pct) && ch > ch_95pct)

    list(n_cells = n, embedding_dim = d, ari_vs_random = ari, ch_index = ch,
         ch_noise_95pct = as.numeric(ch_95pct),
         phantom_flag_ari = flag_ari, phantom_flag_ch_literal = flag_ch_lit,
         phantom_flag_ch_relative = flag_ch_rel,
         phantom_flag_any_relative = flag_ari || flag_ch_rel,
         status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(n_cells = NA, embedding_dim = NA, ari_vs_random = NA, ch_index = NA, ch_noise_95pct = NA,
         phantom_flag_ari = NA, phantom_flag_ch_literal = NA, phantom_flag_ch_relative = NA,
         phantom_flag_any_relative = NA, status = "error", error_msg = conditionMessage(e))
  })

  out[names(res)] <- res
  out
}

results_list <- lapply(seq_len(nrow(nc_manifest)), process_row)
results_df <- do.call(rbind.data.frame, results_list)

n_errors <- sum(results_df$status == "error")
cat(sprintf("\nRows processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) print(table(results_df$error_msg))

cat("\n=== Literal spec threshold (CH > 10) ===\n")
cat("Flagged:", sum(results_df$phantom_flag_ch_literal, na.rm=TRUE), "/", nrow(results_df), "\n")

cat("\n=== Empirically-calibrated relative threshold (CH > noise 95th pct) ===\n")
cat("Flagged:", sum(results_df$phantom_flag_ch_relative, na.rm=TRUE), "/", nrow(results_df), "\n")
cat("\nBy method:\n")
print(aggregate(phantom_flag_ch_relative ~ method, data = results_df, FUN = mean))

cat("\n=== EDA Checkpoint 3 hypothesis: SCTransform v2 / GLM-PCA ===\n")
hyp <- results_df[results_df$method %in% c("pca_sctransform_v2", "pca_glmpca"), ]
print(aggregate(cbind(phantom_flag_ch_relative, ch_index) ~ method, data = hyp, FUN = mean))

cat("\n=== By sparsity level (relative flag) ===\n")
print(aggregate(phantom_flag_ch_relative ~ sparsity + method, data = results_df, FUN = mean))

write.csv(results_df, "data/processed/step6_3_phantom_clustering.csv", row.names = FALSE)
saveRDS(noise_baseline, "data/processed/step6_3_noise_baseline_calibration.rds")
cat("\nWritten results and noise-baseline calibration to disk.\n")
