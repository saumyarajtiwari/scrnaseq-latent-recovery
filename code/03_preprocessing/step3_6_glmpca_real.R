suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(glmpca)
})

MAX_GENES <- 3000
PCA_NV    <- 30
SEED      <- 42
REAL_LR   <- 0.0008
OUT_ROOT  <- "data/processed/pca_glmpca"
LOG_FILE  <- file.path(OUT_ROOT, "step3_6_progress.csv")

log_row <- function(row) {
  df <- as.data.frame(row, stringsAsFactors = FALSE)
  write.table(df, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

process_one <- function(sce_path, out_dir, lr) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_glmpca.rds"))
  if (file.exists(out_file)) {
    return(list(timestamp = as.character(Sys.time()), file_path = sce_path,
                out_path = out_file, status = "skipped_existing",
                nv_used = NA, runtime_sec = NA, msg = ""))
  }
  t0 <- Sys.time()
  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)
    if (is.null(rownames(counts_mat))) rownames(counts_mat) <- paste0("gene", seq_len(nrow(counts_mat)))
    if (is.null(colnames(counts_mat))) colnames(counts_mat) <- paste0("cell", seq_len(ncol(counts_mat)))

    n_genes_orig <- nrow(counts_mat)
    gene_totals <- Matrix::rowSums(counts_mat)
    n_zero_genes <- sum(gene_totals == 0)
    counts_mat <- counts_mat[gene_totals > 0, , drop = FALSE]

    hvg_capped <- FALSE
    if (nrow(counts_mat) > MAX_GENES) {
      gene_means <- Matrix::rowMeans(counts_mat)
      top_genes <- order(gene_means, decreasing = TRUE)[1:MAX_GENES]
      counts_mat <- counts_mat[top_genes, , drop = FALSE]
      hvg_capped <- TRUE
    }

    n_cells <- ncol(counts_mat); n_genes <- nrow(counts_mat)
    nv <- min(PCA_NV, n_cells - 1, n_genes - 1)

    set.seed(SEED)
    fit <- glmpca(as.matrix(counts_mat), L = nv, fam = "nb",
                  ctl = list(verbose = FALSE, maxIter = 25, minIter = 2, lr = lr))

    embedding <- as.matrix(fit$factors)
    loadings <- as.matrix(fit$loadings)
    rownames(embedding) <- colnames(counts_mat)
    rownames(loadings) <- rownames(counts_mat)

    out <- list(
      embedding = embedding, loadings = loadings, nv_used = nv, seed = SEED,
      method = "glmpca_nb", max_iter = 25, lr_used = lr,
      family = "negative_binomial", link = "log",
      n_genes_original = n_genes_orig, n_zero_count_genes_removed = n_zero_genes,
      n_genes_used = n_genes, hvg_capped = hvg_capped, source_file = sce_path,
      true_group = sce$true_group, batch_id = sce$batch_id,
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    list(timestamp = as.character(Sys.time()), file_path = sce_path, out_path = out_file,
         status = "done", nv_used = nv, runtime_sec = out$runtime_sec, msg = "")
  }, error = function(e) {
    list(timestamp = as.character(Sys.time()), file_path = sce_path, out_path = out_file,
         status = "error", nv_used = NA, runtime_sec = NA, msg = conditionMessage(e))
  })
  res
}

excluded_real <- readLines("data/real/glmpca_excluded_real_files.csv")
all_real_files <- list(
  baron = "data/real/baron/baron_harmonized.rds",
  muraro = "data/real/muraro/muraro_harmonized.rds",
  pbmc68k = "data/real/pbmc68k/pbmc68k_harmonized.rds",
  segerstolpe = "data/real/segerstolpe/segerstolpe_harmonized.rds",
  tabula_sapiens_lung = "data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds",
  tasic2018 = "data/real/tasic2018/tasic2018_harmonized.rds"
)
real_files <- all_real_files[!names(all_real_files) %in% excluded_real]

for (name in names(real_files)) {
  row <- process_one(real_files[[name]], file.path(OUT_ROOT, "real", name), lr = REAL_LR)
  log_row(row)
  cat("Real:", name, "-", row$status, "\n")
}
cat("\n=== REAL-DATA PHASE COMPLETE ===\n")
