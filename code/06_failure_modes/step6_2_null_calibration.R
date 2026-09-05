# Item 6: null-calibration for Step 6.2 (Cluster Collapse), silhouette<0.2.
# Unlike Step 6.1's AMI (label-agreement statistic, calibrated via label
# permutation), silhouette is geometry-dependent -- calibrated here via
# Step 6.3's noise-generation template: random Gaussian embeddings at the
# real (n, k, dim) structure, roughly-equal group sizes, same
# simplified_silhouette() function as step4_4to6_secondary_metrics.R.

set.seed(42)
N_DRAWS <- 200

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

d <- read.csv("data/processed/step6_2_cluster_collapse.csv", stringsAsFactors = FALSE)
manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
d <- merge(d, manifest[, c("data_type","source","run_id","method","nv_used")],
           by = c("data_type","source","run_id","method"))

buckets <- unique(d[!is.na(d$n_labeled) & !is.na(d$n_groups_used) & !is.na(d$nv_used),
                     c("n_labeled","n_groups_used","nv_used")])
cat("Calibration buckets:", nrow(buckets), "\n")

calibrate_sil_null <- function(n, k, dim, n_draws) {
  set.seed(2000 + n + k * 7 + dim * 13)
  replicate(n_draws, {
    emb <- matrix(rnorm(n * dim), nrow = n)
    labels <- sample(rep(1:k, length.out = n))
    simplified_silhouette(emb, labels)
  })
}

calib <- list()
for (i in seq_len(nrow(buckets))) {
  n <- buckets$n_labeled[i]; k <- buckets$n_groups_used[i]; dim <- buckets$nv_used[i]
  sil_null <- calibrate_sil_null(n, k, dim, N_DRAWS)
  calib[[paste(n, k, dim)]] <- list(
    n_labeled = n, n_groups_used = k, nv_used = dim,
    null_mean = mean(sil_null), null_median = median(sil_null),
    null_5pct = quantile(sil_null, 0.05), null_1pct = quantile(sil_null, 0.01),
    null_min = min(sil_null)
  )
  if (i %% 5 == 0) cat("Progress:", i, "/", nrow(buckets), "\n")
}
calib_df <- do.call(rbind.data.frame, calib)
cat("\nSilhouette null calibration summary:\n"); print(calib_df)

d$key <- paste(d$n_labeled, d$n_groups_used, d$nv_used)
lookup_5pct <- setNames(calib_df$null_5pct, paste(calib_df$n_labeled, calib_df$n_groups_used, calib_df$nv_used))
lookup_median <- setNames(calib_df$null_median, paste(calib_df$n_labeled, calib_df$n_groups_used, calib_df$nv_used))
d$sil_null_5pct <- lookup_5pct[d$key]
d$sil_null_median <- lookup_median[d$key]

# Relative flag: is the real silhouette below what pure noise typically produces?
# (below the null's 5th percentile = genuinely, surprisingly worse than noise)
d$silhouette_trigger_relative <- d$silhouette < d$sil_null_5pct
d$key <- NULL

write.csv(d, "data/processed/step6_2_cluster_collapse_calibrated.csv", row.names = FALSE)
write.csv(calib_df, "data/processed/step6_2_silhouette_null_calibration.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Literal flag (silhouette<0.2) count:", sum(d$silhouette_trigger, na.rm=TRUE), "\n")
cat("Relative flag (below null 5th pct) count:", sum(d$silhouette_trigger_relative, na.rm=TRUE), "\n")
cat("\nNull median silhouette range across buckets:", range(calib_df$null_median), "\n")
cat("Real data median silhouette (from Step 6.2):", median(d$silhouette, na.rm=TRUE), "\n")
