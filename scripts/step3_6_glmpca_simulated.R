suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(glmpca)
  library(parallel)
})

MAX_GENES <- 3000
PCA_NV <- 30
SEED <- 42
MIN_FREE_GB <- 1.0
BATCH_SIZE <- 100
N_WORKERS <- 6
OUT_ROOT <- "data/processed/pca_glmpca"
LOG_FILE <- file.path(OUT_ROOT, "step3_6_progress.csv")
DISK_MOUNT <- "/mnt/extra2"

log_rows <- function(rows_list) {
  rows_list <- rows_list[!sapply(rows_list, is.null)]
  if (length(rows_list) == 0) return(invisible())
  df <- do.call(rbind, lapply(rows_list, as.data.frame, stringsAsFactors = FALSE))
  write.table(df, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

get_free_gb <- function() {
  out <- system(paste("df --output=avail", DISK_MOUNT, "| tail -1"), intern = TRUE)
  as.numeric(trimws(out)) / 1e6
}

process_one <- function(sce_path, out_dir) {
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
                  ctl = list(verbose = FALSE, maxIter = 25, minIter = 2))

    embedding <- as.matrix(fit$factors)
    loadings <- as.matrix(fit$loadings)
    rownames(embedding) <- colnames(counts_mat)
    rownames(loadings) <- rownames(counts_mat)

    out <- list(
      embedding = embedding, loadings = loadings, nv_used = nv, seed = SEED,
      method = "glmpca_nb", max_iter = 25, lr_used = 0.1,
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

excluded_splatter <- readLines("data/simulated/splatter_excluded_run_ids.csv")
simulator_dirs <- list.dirs("data/simulated/sce", recursive = FALSE)
sim_files <- list()
for (d in simulator_dirs) {
  sim_name <- basename(d)
  files <- list.files(d, full.names = TRUE, pattern = "\\.rds$")
  if (sim_name == "splatter") {
    files <- files[!basename(files) %in% excluded_splatter]
  }
  sim_files[[sim_name]] <- files
}

total_files <- sum(sapply(sim_files, length))
cat("Total simulated files to process:", total_files, "\n\n")

processed_count <- 0
for (sim_name in names(sim_files)) {
  files <- sim_files[[sim_name]]
  out_dir <- file.path(OUT_ROOT, "simulated", sim_name)
  cat("=== Processing", sim_name, "(", length(files), "files ) ===\n")

  n_batches <- ceiling(length(files) / BATCH_SIZE)
  for (b in seq_len(n_batches)) {
    idx_start <- (b - 1) * BATCH_SIZE + 1
    idx_end <- min(b * BATCH_SIZE, length(files))
    batch_files <- files[idx_start:idx_end]

    batch_results <- mclapply(batch_files, function(f) process_one(f, out_dir),
                               mc.cores = N_WORKERS, mc.preschedule = FALSE)
    log_rows(batch_results)

    processed_count <- processed_count + length(batch_files)
    fg <- get_free_gb()
    cat("  [progress] processed:", processed_count, "/", total_files,
        " free disk (GB):", round(fg, 2), "\n")
    if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free. Aborting.")
  }
}
cat("\n=== SIMULATED-DATA PHASE COMPLETE ===\n")
