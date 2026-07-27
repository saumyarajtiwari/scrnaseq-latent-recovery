suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(sctransform)
  library(irlba)
  library(parallel)
})

MAX_GENES    <- 3000
PCA_NV       <- 30
SEED         <- 42
MIN_FREE_GB  <- 1.0
BATCH_SIZE   <- 200
N_WORKERS    <- 2
OUT_ROOT     <- "data/processed/pca_sctransform"
LOG_FILE     <- file.path(OUT_ROOT, "step3_5_progress.csv")
DISK_MOUNT   <- "/mnt/extra2"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

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

process_one <- function(sce_path, out_dir, flags = list()) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_sctpca.rds"))
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
    hvg_capped <- FALSE
    if (n_genes_orig > MAX_GENES) {
      gene_means <- Matrix::rowMeans(counts_mat)
      top_genes <- order(gene_means, decreasing = TRUE)[1:MAX_GENES]
      counts_mat <- counts_mat[top_genes, ]
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
      embedding = embedding,
      singular_values = pca_fit$d,
      var_explained = var_explained,
      nv_used = nv,
      seed = SEED,
      method = "sctransform_v2_pearson_pca",
      latent_var = "log_umi",
      res_clip_range = c(-sqrt(n_cells), sqrt(n_cells)),
      n_genes_original = n_genes_orig,
      n_genes_used = n_genes,
      hvg_capped = hvg_capped,
      centered = TRUE,
      source_file = sce_path,
      true_group = sce$true_group,
      batch_id = sce$batch_id,
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      flags = flags,
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
cat("Splatter exclusion manifest loaded:", length(excluded_splatter), "files\n")

simulator_dirs <- list.dirs("data/simulated/sce", recursive = FALSE)
sim_files <- list()
for (d in simulator_dirs) {
  sim_name <- basename(d)
  files <- list.files(d, full.names = TRUE, pattern = "\\.rds$")
  if (sim_name == "splatter") {
    before <- length(files)
    files <- files[!basename(files) %in% excluded_splatter]
    cat("Splatter: excluded", before - length(files), "files per manifest\n")
  }
  sim_files[[sim_name]] <- files
}

real_files <- list(
  baron         = "data/real/baron/baron_harmonized.rds",
  muraro        = "data/real/muraro/muraro_harmonized.rds",
  pbmc68k       = "data/real/pbmc68k/pbmc68k_harmonized.rds",
  segerstolpe   = "data/real/segerstolpe/segerstolpe_harmonized.rds",
  tabula_sapiens_lung = "data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds",
  tasic2018     = "data/real/tasic2018/tasic2018_harmonized.rds"
)

total_files <- sum(sapply(sim_files, length)) + length(real_files)
cat("Total files to process:", total_files, "\n\n")

## ---- Real datasets: sequential ----
for (name in names(real_files)) {
  row <- process_one(real_files[[name]], file.path(OUT_ROOT, "real", name))
  log_rows(list(row))
  cat("Real:", name, "-", row$status, "\n")
}

## ---- Simulated datasets: parallel batches, matching Step 1.3's validated config ----
processed_count <- length(real_files)
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

cat("\n=== STEP 3.5 FULL RUN COMPLETE ===\n")
