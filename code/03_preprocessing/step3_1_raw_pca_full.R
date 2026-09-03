suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

RAW_PCA_NV   <- 50
SEED         <- 42
MIN_FREE_GB  <- 1.0
CHECK_DISK_EVERY <- 500
LOG_FILE     <- "data/processed/pca_raw/step3_1_progress.csv"

dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)

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
  out <- system("df --output=avail /mnt/extra | tail -1", intern = TRUE)
  as.numeric(trimws(out)) / 1e6
}

run_raw_pca_one <- function(sce_path, out_dir, flags = list()) {
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(sce_path)), "_rawpca.rds"))
  if (file.exists(out_file)) return(list(status = "skipped_existing"))

  t0 <- Sys.time()
  res <- tryCatch({
    sce <- readRDS(sce_path)
    counts_mat <- counts(sce)

    n_cells <- ncol(counts_mat); n_genes <- nrow(counts_mat)
    nv <- min(RAW_PCA_NV, n_cells - 1, n_genes - 1)

    set.seed(SEED)
    pca_fit <- irlba(t(counts_mat), nv = nv, center = FALSE, scale = FALSE)
    embedding <- pca_fit$u %*% diag(pca_fit$d)
    rownames(embedding) <- colnames(counts_mat)
    var_explained <- (pca_fit$d^2) / sum(pca_fit$d^2)

    out <- list(
      embedding = embedding,
      singular_values = pca_fit$d,
      var_explained = var_explained,
      nv_used = nv,
      seed = SEED,
      method = "irlba_raw_pca",
      source_file = sce_path,
      true_group = sce$true_group,
      batch_id = sce$batch_id,
      sim_params = tryCatch(metadata(sce)$run_params, error = function(e) NULL),
      flags = flags,
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    rm(sce, counts_mat, pca_fit, embedding); gc(verbose = FALSE)
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

manifest_candidates <- c(
  "data/simulated/splatter_excluded_run_ids.csv",
  "data/simulated/excluded_files.csv",
  "data/simulated/sce/splatter/excluded_run_ids.csv"
)
excluded_splatter <- character(0)
manifest_found <- FALSE
for (m in manifest_candidates) {
  if (file.exists(m)) {
    excluded_splatter <- readLines(m)
    manifest_found <- TRUE
    cat("Exclusion manifest found:", m, "-", length(excluded_splatter), "files excluded\n")
    break
  }
}
if (!manifest_found) {
  cat("*** WARNING: no Splatter exclusion manifest found at any candidate path.\n")
  cat("*** Proceeding WITHOUT excluding the 6 known unresolved files.\n")
  cat("*** This is logged, not silent -- flag to Claude if you have the manifest path.\n\n")
}

simulator_dirs <- list.dirs("data/simulated/sce", recursive = FALSE)
sim_files <- list()
for (d in simulator_dirs) {
  sim_name <- basename(d)
  files <- list.files(d, full.names = TRUE, pattern = "\\.rds$")
  if (sim_name == "splatter" && length(excluded_splatter) > 0) {
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
cat("Total files to process (excluding already-done on resume):", total_files, "\n\n")

counter <- 0
for (name in names(real_files)) {
  counter <- counter + 1
  run_raw_pca_one(real_files[[name]], file.path("data/processed/pca_raw/real", name))
  if (counter %% CHECK_DISK_EVERY == 0) {
    fg <- get_free_gb()
    if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free. Aborting.")
  }
}

for (sim_name in names(sim_files)) {
  files <- sim_files[[sim_name]]
  out_dir <- file.path("data/processed/pca_raw/simulated", sim_name)
  cat("=== Processing", sim_name, "(", length(files), "files) ===\n")

  for (i in seq_along(files)) {
    counter <- counter + 1
    flags <- list()
    if (sim_name == "scdesign3") flags$note <- "depth500 underdelivery not excluded, flagged if applicable"
    run_raw_pca_one(files[i], out_dir, flags = flags)

    if (counter %% CHECK_DISK_EVERY == 0) {
      fg <- get_free_gb()
      cat("  [progress] processed:", counter, "/", total_files, " free disk (GB):", round(fg, 2), "\n")
      if (fg < MIN_FREE_GB) stop("Disk guard triggered: only ", fg, "GB free. Aborting at file ", i, " of ", sim_name)
    }
  }
}

cat("\n=== STEP 3.1 FULL RUN COMPLETE ===\n")
