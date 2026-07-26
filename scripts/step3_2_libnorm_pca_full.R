suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

TARGET_SUM   <- 10000
PCA_NV       <- 30
SEED         <- 42
MIN_FREE_GB  <- 1.0
CHECK_DISK_EVERY <- 500
OUT_ROOT     <- "data/processed/pca_libnorm"
LOG_FILE     <- file.path(OUT_ROOT, "step3_2_progress.csv")
DISK_MOUNT   <- "/mnt/extra2"

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

log_row <- function(file_path, out_path, status, nv_used, runtime_sec, msg = "") {
  row <- data.frame(
    timestamp   = as.character(Sys.time()),
    file_path   = file_path,
    out_path    = out_path,
    status      = status,
    nv_used     = nv_used,
    runtime_sec = runtime_sec,
    msg         = msg,
    stringsAsFactors = FALSE
  )
  write.table(row, LOG_FILE, sep = ",", row.names = FALSE,
              col.names = !file.exists(LOG_FILE), append = file.exists(LOG_FILE))
}

get_free_gb <- function() {
  out <- system(paste("df --output=avail", DISK_MOUNT, "| tail -1"), intern = TRUE)
  as.numeric(trimws(out)) / 1e6
}

run_libnorm_pca_one <- function(sce_path, out_dir, flags = list()) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_libnormpca.rds"))
  if (file.exists(out_file)) return(list(status = "skipped_existing"))

  t0 <- Sys.time()
  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)

    lib_sizes <- Matrix::colSums(counts_mat)
    n_zero <- sum(lib_sizes == 0)
    safe_lib_sizes <- lib_sizes
    safe_lib_sizes[safe_lib_sizes == 0] <- 1  # 0-count cells stay 0 after division; avoids div-by-zero

    scale_factor <- TARGET_SUM / safe_lib_sizes
    normalized <- counts_mat %*% Diagonal(x = scale_factor)  # genes x cells, sparsity-preserving

    n_cells <- ncol(normalized); n_genes <- nrow(normalized)
    nv <- min(PCA_NV, n_cells - 1, n_genes - 1)

    set.seed(SEED)
    pca_fit <- irlba(t(normalized), nv = nv, center = TRUE, scale = FALSE)
    embedding <- pca_fit$u %*% diag(pca_fit$d)
    rownames(embedding) <- colnames(counts_mat)
    var_explained <- (pca_fit$d^2) / sum(pca_fit$d^2)

    out <- list(
      embedding = embedding,
      singular_values = pca_fit$d,
      var_explained = var_explained,
      nv_used = nv,
      seed = SEED,
      method = "irlba_libnorm_pca",
      target_sum = TARGET_SUM,
      centered = TRUE,
      n_zero_libsize_cells = n_zero,
      source_file = sce_path,
      true_group = sce$true_group,
      batch_id = sce$batch_id,
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      flags = flags,
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    rm(sce, counts_mat, normalized, pca_fit, embedding); gc(verbose = FALSE)
    list(status = "done", nv_used = nv, runtime = out$runtime_sec)
  }, error = function(e) {
    list(status = "error", msg = conditionMessage(e))
  })

  if (res$status == "done") {
    log_row(sce_path, out_file, "done", res$nv_used, res$runtime)
  } else if (res$status == "error") {
    log_row(sce_path, out_file, "error", NA, NA, res$msg)
    cat("  ERROR:", sce_path, "-", res$msg, "\n")
  }
  res
}

## ---- Exclusion manifest (now resolved) ----
excluded_splatter <- readLines("data/simulated/splatter_excluded_run_ids.csv")
cat("Splatter exclusion manifest loaded:", length(excluded_splatter), "files\n")

## ---- Build file list ----
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

counter <- 0
for (name in names(real_files)) {
  counter <- counter + 1
  run_libnorm_pca_one(real_files[[name]], file.path(OUT_ROOT, "real", name))
  if (counter %% CHECK_DISK_EVERY == 0) {
    fg <- get_free_gb()
    if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free. Aborting.")
  }
}

for (sim_name in names(sim_files)) {
  files <- sim_files[[sim_name]]
  out_dir <- file.path(OUT_ROOT, "simulated", sim_name)
  cat("=== Processing", sim_name, "(", length(files), "files) ===\n")

  for (i in seq_along(files)) {
    counter <- counter + 1
    run_libnorm_pca_one(files[i], out_dir)

    if (counter %% CHECK_DISK_EVERY == 0) {
      fg <- get_free_gb()
      cat("  [progress] processed:", counter, "/", total_files, " free disk (GB):", round(fg, 2), "\n")
      if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free. Aborting at file ", i, " of ", sim_name)
    }
  }
}

cat("\n=== STEP 3.2 FULL RUN COMPLETE ===\n")
