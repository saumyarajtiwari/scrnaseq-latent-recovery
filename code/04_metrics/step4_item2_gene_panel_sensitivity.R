suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

MAX_GENES <- 3000
SEED <- 42
TARGET_SUM <- 10000

manifest <- read.csv("data/processed/step3_10_pre_internal_method_backup/embedding_manifest.csv", stringsAsFactors = FALSE)
sub <- manifest[manifest$data_type == "simulated" & manifest$method == "pca_raw" &
                 !(manifest$is_null_control %in% TRUE) &
                 manifest$dropout == "low" & manifest$batch == "none" &
                 manifest$clipping == "none" & manifest$n_cells == 1000 &
                 manifest$gene_strategy == "all", ]
groups <- unique(sub[, c("source", "run_id")])
cat("Total (simulator, run_id) groups to process:", nrow(groups), "\n")

principal_angles_deg <- function(true_basis, est_basis) {
  Qt <- qr.Q(qr(true_basis)); Qe <- qr.Q(qr(est_basis))
  s <- svd(t(Qt) %*% Qe)$d
  s <- pmin(pmax(s, -1), 1)
  acos(s) * 180 / pi
}

apply_hvg_cap <- function(counts_mat) {
  if (nrow(counts_mat) <= MAX_GENES) return(counts_mat)
  gene_means <- Matrix::rowMeans(counts_mat)
  top_genes <- order(gene_means, decreasing = TRUE)[1:MAX_GENES]
  subset_check <- counts_mat[top_genes, , drop = FALSE]
  zero_after_subset <- which(Matrix::colSums(subset_check) == 0)
  if (length(zero_after_subset) > 0) {
    rescue_genes <- sapply(zero_after_subset, function(ci) which.max(counts_mat[, ci]))
    top_genes <- union(top_genes, unique(rescue_genes))
  }
  counts_mat[top_genes, , drop = FALSE]
}

fit_method <- function(counts_mat, method) {
  if (method == "pca_raw") {
    nv <- min(50, ncol(counts_mat) - 1, nrow(counts_mat) - 1)
    set.seed(SEED)
    fit <- irlba(t(counts_mat), nv = nv, center = FALSE, scale = FALSE)
    return(list(loadings = fit$v))
  }
  lib_sizes <- Matrix::colSums(counts_mat)
  safe_lib <- lib_sizes; safe_lib[safe_lib == 0] <- 1
  normalized <- counts_mat %*% Diagonal(x = TARGET_SUM / safe_lib)

  if (method == "pca_libnorm") {
    mat <- normalized
  } else if (method == "pca_log") {
    mat <- log1p(normalized)
  } else if (method == "pca_shiftedlog") {
    delta_opt <- median(normalized@x)
    mat <- normalized
    mat@x <- log(mat@x + delta_opt) - log(delta_opt)
  }
  nv <- min(30, ncol(mat) - 1, nrow(mat) - 1)
  set.seed(SEED)
  fit <- irlba(t(mat), nv = nv, center = TRUE, scale = FALSE)
  list(loadings = fit$v)
}

METHODS <- c("pca_raw", "pca_libnorm", "pca_log", "pca_shiftedlog")
results <- list()

for (i in seq_len(nrow(groups))) {
  sim_name <- groups$source[i]; rid <- groups$run_id[i]
  sce_path <- sprintf("data/simulated/sce/%s/%s_sce_run_%05d.rds", sim_name, sim_name, rid)

  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)
    tgm_full <- S4Vectors::metadata(sce)$true_group_means
    K <- ncol(tgm_full); r <- K - 1

    counts_capped <- apply_hvg_cap(counts_mat)
    hvg_genes <- rownames(counts_capped)

    tgm_sub <- tgm_full[rownames(tgm_full) %in% hvg_genes, , drop = FALSE]
    tgm_sub <- tgm_sub[match(hvg_genes, rownames(tgm_sub)), , drop = FALSE]
    tgm_centered <- tgm_sub - rowMeans(tgm_sub)
    true_basis <- svd(tgm_centered)$u[, 1:r, drop = FALSE]

    for (method in METHODS) {
      m_res <- tryCatch({
        fit <- fit_method(counts_capped, method)
        est_basis <- fit$loadings[, 1:r, drop = FALSE]
        angles <- principal_angles_deg(true_basis, est_basis)
        angles_rad <- angles * pi / 180
        list(status = "ok",
             subspace_recovery_score_samegenes = mean(cos(angles_rad)^2),
             grassmann_distance_samegenes = sqrt(sum(sin(angles_rad)^2)),
             n_genes_used_samegenes = nrow(counts_capped))
      }, error = function(e) list(status = "error",
           subspace_recovery_score_samegenes = NA, grassmann_distance_samegenes = NA,
           n_genes_used_samegenes = NA))

      results[[length(results) + 1]] <- data.frame(
        source = sim_name, run_id = rid, method = method,
        status = m_res$status,
        subspace_recovery_score_samegenes = m_res$subspace_recovery_score_samegenes,
        grassmann_distance_samegenes = m_res$grassmann_distance_samegenes,
        n_genes_used_samegenes = m_res$n_genes_used_samegenes,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    for (method in METHODS) {
      results[[length(results) + 1]] <<- data.frame(
        source = sim_name, run_id = rid, method = method, status = "group_error",
        subspace_recovery_score_samegenes = NA, grassmann_distance_samegenes = NA,
        n_genes_used_samegenes = NA, stringsAsFactors = FALSE)
    }
  })

  if (i %% 20 == 0) cat("Progress:", i, "/", nrow(groups), "\n")
}

final <- do.call(rbind, results)

# Join against existing full-gene scores for direct comparison
full_scores <- read.csv("results/step4_metrics/step4_1to3_subspace_metrics.csv", stringsAsFactors = FALSE)
final <- merge(final, full_scores[, c("source","run_id","method","subspace_recovery_score","grassmann_distance")],
               by = c("source","run_id","method"), all.x = TRUE)
final$recovery_score_delta <- final$subspace_recovery_score_samegenes - final$subspace_recovery_score

write.csv(final, "data/processed/step4_item2_gene_panel_sensitivity.csv", row.names = FALSE)
cat("\n=== ITEM 2 SENSITIVITY CHECK COMPLETE ===\n")
cat("Total rows:", nrow(final), "\n")
cat("\nMean recovery_score_delta (samegenes - full) by method:\n")
print(aggregate(recovery_score_delta ~ method, data = final, FUN = mean, na.rm = TRUE))
