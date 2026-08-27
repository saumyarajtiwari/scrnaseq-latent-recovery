# Step 6.3 — Detect Phantom Clustering
# k-means (k=3) on null-control embeddings (single true population, n_groups=1).
# Two tests, reported separately given the interpretation gap flagged before coding:
#   (a) ARI vs. a single random label draw, exactly as specified. Note: ARI is
#       chance-corrected by construction, so its expected value against any
#       independent random partition is ~0 regardless of whether the k-means
#       clustering reflects real structure; a single draw crossing 0.2 is
#       largely a function of sampling luck at small-to-moderate n, not a
#       reliable structure/no-structure discriminant. Reported faithfully,
#       not treated as the primary indicator.
#   (b) Calinski-Harabasz index (>10 threshold, per spec) as the primary,
#       non-chance-corrected indicator of genuine artifactual cluster structure
#       in a known single-population dataset.
# Directly tests the EDA Checkpoint 3 pre-registered hypothesis: SCTransform v2
# and GLM-PCA's null-control panels showed visually non-unimodal structure,
# flagged there as "not yet confirmed as a real finding."

suppressPackageStartupMessages(library(aricode))

set.seed(42)

nc_manifest <- read.csv("data/processed/step6_3_null_control_manifest.csv", stringsAsFactors = FALSE)
cat("Total null-control files:", nrow(nc_manifest), "\n")

calinski_harabasz <- function(emb, cluster) {
  k <- length(unique(cluster))
  n <- nrow(emb)
  overall_centroid <- colMeans(emb)
  between_ss <- 0
  within_ss <- 0
  for (c in unique(cluster)) {
    pts <- emb[cluster == c, , drop = FALSE]
    n_c <- nrow(pts)
    centroid_c <- colMeans(pts)
    between_ss <- between_ss + n_c * sum((centroid_c - overall_centroid)^2)
    within_ss <- within_ss + sum(sweep(pts, 2, centroid_c, "-")^2)
  }
  if (within_ss == 0) return(Inf)
  (between_ss / (k - 1)) / (within_ss / (n - k))
}

process_row <- function(i) {
  r <- nc_manifest[i, ]
  out <- list(
    method = r$method, source = r$source, run_id = r$run_id, replicate = r$replicate,
    sparsity = r$sparsity, n_cells = NA, ari_vs_random = NA, ch_index = NA,
    phantom_flag_ari = NA, phantom_flag_ch = NA, phantom_flag_any = NA,
    status = "ok", error_msg = NA_character_
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    n <- nrow(emb)

    set.seed(42)
    km <- kmeans(emb, centers = 3, nstart = 25)$cluster

    set.seed(42 + r$run_id * 10 + r$replicate)  # reproducible but distinct per file
    random_labels <- sample(1:3, n, replace = TRUE)
    ari <- aricode::ARI(as.integer(km), as.integer(random_labels))

    ch <- calinski_harabasz(emb, km)

    list(n_cells = n, ari_vs_random = ari, ch_index = ch,
         phantom_flag_ari = isTRUE(ari > 0.2), phantom_flag_ch = isTRUE(ch > 10),
         phantom_flag_any = isTRUE(ari > 0.2) || isTRUE(ch > 10),
         status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(n_cells = NA, ari_vs_random = NA, ch_index = NA,
         phantom_flag_ari = NA, phantom_flag_ch = NA, phantom_flag_any = NA,
         status = "error", error_msg = conditionMessage(e))
  })

  out[names(res)] <- res
  out
}

results_list <- lapply(seq_len(nrow(nc_manifest)), process_row)
results_df <- do.call(rbind.data.frame, results_list)

n_errors <- sum(results_df$status == "error")
cat(sprintf("Rows processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) print(table(results_df$error_msg))

cat("\n=== ARI-vs-random-label test ===\n")
print(summary(results_df$ari_vs_random))
cat("Flagged (ARI > 0.2):", sum(results_df$phantom_flag_ari, na.rm=TRUE), "/", nrow(results_df), "\n")

cat("\n=== Calinski-Harabasz test ===\n")
print(summary(results_df$ch_index))
cat("Flagged (CH > 10):", sum(results_df$phantom_flag_ch, na.rm=TRUE), "/", nrow(results_df), "\n")

cat("\n=== By method: CH-flagged rate (the primary indicator) ===\n")
print(aggregate(phantom_flag_ch ~ method, data = results_df, FUN = mean))

cat("\n=== EDA Checkpoint 3 hypothesis check: SCTransform v2 / GLM-PCA specifically ===\n")
flagged_methods <- results_df[results_df$method %in% c("pca_sctransform_v2", "pca_glmpca"), ]
print(aggregate(cbind(phantom_flag_ch, ch_index) ~ method, data = flagged_methods, FUN = mean))

write.csv(results_df, "data/processed/step6_3_phantom_clustering.csv", row.names = FALSE)
cat("\nWritten to data/processed/step6_3_phantom_clustering.csv\n")
