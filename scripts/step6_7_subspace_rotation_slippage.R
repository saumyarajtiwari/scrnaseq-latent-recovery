# Step 6.7 — Detect Subspace Rotation Slippage (v2: gene-subset fix)
#
# [Same header as before regarding methodology and validation]
#
# BUG FIX (found in dry run: 167/500 errors, "non-conformable arguments"):
# 82 pca_glmpca + 85 pca_sctransform_v2 = 167 exactly. Both scripts cap
# genes at MAX_GENES=3000 (confirmed in step3_5_sctransform_with_loadings.R
# and step3_6_glmpca_simulated.R) -- for Splatter's 10,000-gene files this
# triggers HVG subsetting, giving loadings with ~3000 genes while
# true_group_means retains the full gene set. NOT a gene_strategy grid-axis
# effect (verified directly: an hvg2000-labeled pca_raw row had full
# 10,000-gene loadings matching true_group_means exactly -- gene_strategy
# is metadata only for the 4 linear methods, no actual subsetting occurs).
# Fixed by subsetting true_group_means to loadings' actual rownames
# (matched by gene name, not position) before computing Q_true. Since
# SCTransform v2 and GLM-PCA each select their own internal HVG subset
# (potentially different from each other), the true-subspace cache can no
# longer be shared across all 6 methods per source file -- now keyed by
# (source_file, method).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args
set.seed(42)

N_PCS_CHECK <- 3
ANGLE_THRESHOLD <- 30

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
work <- manifest[manifest$data_type == "simulated" &
                  !(manifest$is_null_control %in% TRUE) & manifest$status == "ok", ]
cat(sprintf("Total candidate rows: %d\n", nrow(work)))

if (dry_run) {
  set.seed(1)
  work <- work[sample(nrow(work), 500), ]
  cat(sprintf("DRY RUN: subset to %d rows\n", nrow(work)))
}

true_subspace_cache_dir <- "data/processed/step6_true_subspace_cache"
dir.create(true_subspace_cache_dir, showWarnings = FALSE, recursive = TRUE)

get_true_subspace <- function(source_file, gene_names, method) {
  cache_key <- paste0(gsub("[/\\\\]", "_", source_file), "__", method)
  cache_path <- file.path(true_subspace_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(cache_path)) {
    result <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(result)) return(result)
  }

  src <- readRDS(source_file)
  tgm <- S4Vectors::metadata(src)$true_group_means

  # Subset to the SAME genes actually used in this file's loadings,
  # matched by name -- required for SCTransform v2/GLM-PCA's internal
  # HVG capping, harmless no-op for the 4 methods that use all genes.
  common_genes <- intersect(gene_names, rownames(tgm))
  if (length(common_genes) < length(gene_names)) {
    stop(sprintf("only %d/%d loadings genes found in true_group_means",
                  length(common_genes), length(gene_names)))
  }
  tgm <- tgm[gene_names, , drop = FALSE]  # match loadings' exact gene order

  tgm_centered <- sweep(tgm, 1, rowMeans(tgm), "-")
  sv <- svd(tgm_centered)
  tol <- max(dim(tgm_centered)) * .Machine$double.eps * max(sv$d)
  rank_true <- sum(sv$d > tol)
  Q_true <- sv$u[, seq_len(rank_true), drop = FALSE]

  result <- list(Q_true = Q_true, rank_true = rank_true)
  tmp_path <- paste0(cache_path, ".tmp_", Sys.getpid())
  saveRDS(result, tmp_path)
  file.rename(tmp_path, cache_path)
  result
}

process_row <- function(i) {
  r <- work[i, ]
  out <- list(
    source = r$source, run_id = r$run_id, method = r$method, file_path = r$file_path,
    rank_true = NA, gene_strategy = r$gene_strategy, sparsity = r$sparsity, depth = r$depth,
    dropout = r$dropout, separability = r$separability, batch = r$batch, clipping = r$clipping,
    angle_pc1 = NA, angle_pc2 = NA, angle_pc3 = NA,
    flag_pc1 = NA, flag_pc2 = NA, flag_pc3 = NA, first_pc_flagged = NA_integer_,
    status = "ok", error_msg = NA_character_
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    loadings <- obj$loadings
    n_pcs_avail <- min(N_PCS_CHECK, ncol(loadings))
    gene_names <- rownames(loadings)
    if (is.null(gene_names)) stop("loadings has no rownames (gene IDs)")

    ts <- get_true_subspace(obj$source_file, gene_names, r$method)
    Q_true <- ts$Q_true

    angles <- rep(NA_real_, N_PCS_CHECK)
    for (k in seq_len(n_pcs_avail)) {
      v <- loadings[, k]
      v_norm <- sqrt(sum(v^2))
      if (v_norm == 0) next
      v <- v / v_norm
      proj_norm <- sqrt(sum((crossprod(Q_true, v))^2))
      proj_norm <- min(proj_norm, 1)
      angles[k] <- acos(proj_norm) * 180 / pi
    }

    flags <- angles > ANGLE_THRESHOLD
    first_flag <- which(flags)
    first_flag <- if (length(first_flag) > 0) min(first_flag) else NA_integer_

    list(rank_true = ts$rank_true,
         angle_pc1 = angles[1], angle_pc2 = angles[2], angle_pc3 = angles[3],
         flag_pc1 = flags[1], flag_pc2 = flags[2], flag_pc3 = flags[3],
         first_pc_flagged = first_flag, status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(rank_true = NA, angle_pc1 = NA, angle_pc2 = NA, angle_pc3 = NA,
         flag_pc1 = NA, flag_pc2 = NA, flag_pc3 = NA, first_pc_flagged = NA_integer_,
         status = "error", error_msg = conditionMessage(e))
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

cat("\nAngle summaries:\n")
print(summary(results_df[, c("angle_pc1","angle_pc2","angle_pc3")]))
cat("\nFirst-PC-flagged distribution:\n")
print(table(results_df$first_pc_flagged, useNA = "always"))

suffix <- if (dry_run) "_DRYRUN" else ""
write.csv(results_df, paste0("data/processed/step6_7_subspace_rotation_slippage", suffix, ".csv"), row.names = FALSE)
cat("Written.\n")
