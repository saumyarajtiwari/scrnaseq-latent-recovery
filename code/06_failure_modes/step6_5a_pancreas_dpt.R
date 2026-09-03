# Step 6.5a — Over-Smoothing detection, real-data pseudotime branch
# (v2: fixed file-discovery bug that accidentally picked up the stale
# data/processed/pca_sctransform/ directory -- a superseded, pre-loadings-
# backfill output from Step 3.5's earlier rescue attempts, ~3 weeks older
# and half the file size of the current pca_sctransform_v2 output. Only
# the 6 canonical method directories are now used, matching the same
# METHOD_MAP convention as step3_10_manifest.R.)

suppressPackageStartupMessages(library(RSpectra))
set.seed(42)

METHODS <- c("pca_raw", "pca_libnorm", "pca_log", "pca_shiftedlog",
             "pca_sctransform_v2", "pca_glmpca")

lineage_map <- list(
  muraro = list(progenitor = c("duct"), mature = c("alpha","beta","delta","epsilon","pp")),
  baron = list(progenitor = c("ductal"), mature = c("alpha","beta","delta","epsilon","gamma")),
  segerstolpe = list(progenitor = c("ductal cell"),
                      mature = c("alpha cell","beta cell","delta cell","epsilon cell","gamma cell"))
)

diffusion_pseudotime <- function(emb, root_idx, n_eigen = 10) {
  n <- nrow(emb)
  d2 <- as.matrix(dist(emb))^2
  sigma2 <- median(d2[d2 > 0])
  K <- exp(-d2 / sigma2)
  deg <- rowSums(K)
  D_inv_sqrt <- 1 / sqrt(deg)
  K_sym <- sweep(sweep(K, 1, D_inv_sqrt, "*"), 2, D_inv_sqrt, "*")
  K_sym <- (K_sym + t(K_sym)) / 2

  eig <- RSpectra::eigs_sym(K_sym, k = n_eigen + 1, which = "LA")
  lambda <- eig$values[-1]
  phi <- eig$vectors[, -1, drop = FALSE]
  psi <- sweep(phi, 1, D_inv_sqrt, "*")

  lambda <- pmin(pmax(lambda, 1e-8), 1 - 1e-8)
  weights <- lambda / (1 - lambda)

  diffs <- sweep(psi, 2, psi[root_idx, ], "-")
  sqrt(rowSums(sweep(diffs^2, 2, weights^2, "*")))
}

results <- list()

for (ds in names(lineage_map)) {
  cat(sprintf("\n=== %s ===\n", ds))
  for (method in METHODS) {
    fp <- list.files(file.path("data/processed", method, "real", ds), full.names = TRUE)
    if (length(fp) != 1) {
      cat(sprintf("  %s: SKIP (expected 1 file, found %d)\n", method, length(fp)))
      next
    }
    res <- tryCatch({
      obj <- readRDS(fp)
      emb <- obj$embedding
      tg <- as.character(obj$true_group)

      prog <- lineage_map[[ds]]$progenitor
      mat  <- lineage_map[[ds]]$mature
      keep <- tg %in% c(prog, mat)
      n_prog <- sum(tg %in% prog); n_mat <- sum(tg %in% mat)
      if (n_prog < 5 || n_mat < 5) stop(sprintf("insufficient cells: progenitor=%d mature=%d", n_prog, n_mat))

      emb_sub <- emb[keep, , drop = FALSE]
      tg_sub <- tg[keep]
      ref_binary <- as.integer(tg_sub %in% mat)

      prog_idx <- which(tg_sub %in% prog)
      prog_centroid <- colMeans(emb_sub[prog_idx, , drop = FALSE])
      dists <- sqrt(rowSums(sweep(emb_sub[prog_idx, , drop = FALSE], 2, prog_centroid, "-")^2))
      root_idx <- prog_idx[which.min(dists)]

      dpt <- diffusion_pseudotime(emb_sub, root_idx)
      rho <- suppressWarnings(cor(dpt, ref_binary, method = "spearman"))

      list(n_cells_used = length(tg_sub), n_progenitor = n_prog, n_mature = n_mat,
           spearman_dpt_vs_reference = rho, over_smoothing_flag = isTRUE(rho < 0.6),
           status = "ok", error_msg = NA_character_)
    }, error = function(e) {
      list(n_cells_used = NA, n_progenitor = NA, n_mature = NA,
           spearman_dpt_vs_reference = NA, over_smoothing_flag = NA,
           status = "error", error_msg = conditionMessage(e))
    })
    res$source <- ds; res$method <- method; res$file_path <- fp
    results[[length(results) + 1]] <- res
    cat(sprintf("  %s: status=%s rho=%s flag=%s %s\n", method, res$status,
                ifelse(is.na(res$spearman_dpt_vs_reference), "NA", round(res$spearman_dpt_vs_reference, 3)),
                res$over_smoothing_flag, ifelse(res$status == "error", paste("ERROR:", res$error_msg), "")))
  }
}

results_df <- do.call(rbind.data.frame, results)
write.csv(results_df, "data/processed/step6_5a_pancreas_dpt.csv", row.names = FALSE)
cat("\nTotal rows:", nrow(results_df), "(expected 18)\n")
cat("Written to data/processed/step6_5a_pancreas_dpt.csv\n")
