suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(sctransform)
  library(irlba)
})

MAX_GENES <- 3000
PCA_NV    <- 30
SEED      <- 42
OUT_ROOT  <- "data/processed/pca_sctransform"
LOG_FILE  <- file.path(OUT_ROOT, "step3_5_rescue2_progress.csv")

log_row <- function(file_path, out_path, status, nv_used, runtime_sec, msg = "") {
  row <- data.frame(
    timestamp = as.character(Sys.time()), file_path = file_path, out_path = out_path,
    status = status, nv_used = nv_used, runtime_sec = runtime_sec, msg = msg,
    stringsAsFactors = FALSE
  )
  write.table(row, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

process_one <- function(sce_path, out_dir) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_sctpca.rds"))
  if (file.exists(out_file)) return(invisible())

  t0 <- Sys.time()
  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)
    if (is.null(rownames(counts_mat))) rownames(counts_mat) <- paste0("gene", seq_len(nrow(counts_mat)))
    if (is.null(colnames(counts_mat))) colnames(counts_mat) <- paste0("cell", seq_len(ncol(counts_mat)))

    full_lib_sizes <- Matrix::colSums(counts_mat)
    zero_cell_idx <- which(full_lib_sizes == 0)
    excluded_cell_ids <- colnames(counts_mat)[zero_cell_idx]
    excluded_true_group <- sce$true_group[zero_cell_idx]

    if (length(zero_cell_idx) > 0) {
      counts_mat <- counts_mat[, -zero_cell_idx, drop = FALSE]
      kept_true_group <- sce$true_group[-zero_cell_idx]
      kept_batch_id <- sce$batch_id[-zero_cell_idx]
    } else {
      kept_true_group <- sce$true_group
      kept_batch_id <- sce$batch_id
    }

    n_genes_orig <- nrow(counts_mat)
    hvg_capped <- FALSE
    if (n_genes_orig > MAX_GENES) {
      gene_means <- Matrix::rowMeans(counts_mat)
      top_genes <- order(gene_means, decreasing = TRUE)[1:MAX_GENES]

      subset_check <- counts_mat[top_genes, , drop = FALSE]
      zero_after_subset <- which(Matrix::colSums(subset_check) == 0)
      if (length(zero_after_subset) > 0) {
        rescue_genes <- sapply(zero_after_subset, function(ci) which.max(counts_mat[, ci]))
        top_genes <- union(top_genes, unique(rescue_genes))
      }
      counts_mat <- counts_mat[top_genes, , drop = FALSE]
      hvg_capped <- TRUE
    }

    set.seed(SEED)
    vst_out <- sctransform::vst(counts_mat, latent_var = "log_umi", vst.flavor = "v2",
                                 return_gene_attr = FALSE, verbosity = 0)
    residuals <- vst_out$y

    n_cells <- ncol(residuals); n_genes <- nrow(residuals)
    nv <- min(PCA_NV, n_cells - 1, n_genes - 1)

    set.seed(SEED)
    pca_fit <- irlba(t(residuals), nv = nv, center = TRUE, scale = FALSE)
    embedding <- pca_fit$u %*% diag(pca_fit$d)
    rownames(embedding) <- colnames(counts_mat)
    var_explained <- (pca_fit$d^2) / sum(pca_fit$d^2)

    out <- list(
      embedding = embedding, singular_values = pca_fit$d, var_explained = var_explained,
      nv_used = nv, seed = SEED, method = "sctransform_v2_pearson_pca",
      latent_var = "log_umi", res_clip_range = c(-sqrt(n_cells), sqrt(n_cells)),
      n_genes_original = n_genes_orig, n_genes_used = n_genes, hvg_capped = hvg_capped,
      centered = TRUE, source_file = sce_path,
      true_group = kept_true_group, batch_id = kept_batch_id,
      excluded_zero_count_cells = excluded_cell_ids,
      excluded_zero_count_cells_true_group = excluded_true_group,
      exclusion_reason = "SCTransform requires log10(total UMI) as latent variable, undefined for zero-count cells -- genuine methodological limitation, not present in library-size-based methods (3.2-3.4) which degrade gracefully via a zero-substitution guard",
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      flags = list(note = "reprocessed via rescue2 script: native zero-count cell exclusion"),
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    saveRDS(out, out_file)
    list(status = "done", nv_used = nv, runtime = out$runtime_sec, excluded = length(excluded_cell_ids))
  }, error = function(e) list(status = "error", msg = conditionMessage(e)))

  if (res$status == "done") {
    log_row(sce_path, out_file, "done", res$nv_used, res$runtime, paste("excluded_zero_cells:", res$excluded))
    cat(basename(sce_path), "- done, excluded", res$excluded, "zero-count cell(s)\n")
  } else {
    log_row(sce_path, out_file, "error", NA, NA, res$msg)
    cat(basename(sce_path), "- STILL FAILING:", res$msg, "\n")
  }
}

still_failing <- c("splatter_sce_run_01096","splatter_sce_run_01189","splatter_sce_run_01277",
"splatter_sce_run_01546","splatter_sce_run_01907","splatter_sce_run_04382","splatter_sce_run_04428",
"splatter_sce_run_05641","splatter_sce_run_06003","splatter_sce_run_06453","splatter_sce_run_07219",
"splatter_sce_run_07486","splatter_sce_run_07580","splatter_sce_run_07847","splatter_sce_run_09603",
"splatter_sce_run_09646","splatter_sce_run_10817")

for (fn in still_failing) {
  f <- file.path("data/simulated/sce/splatter", paste0(fn, ".rds"))
  process_one(f, "data/processed/pca_sctransform/simulated/splatter")
}

cat("\n=== RESCUE2 RUN COMPLETE ===\n")
