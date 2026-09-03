# =============================================================================
# convert_to_sce.R
#
# Converts production count matrices + cell_meta + run_params + ground truth
# into unified SingleCellExperiment objects -- the actual Step 1.8 deliverable.
#
# Ground truth attachment:
#   Splatter: direct run_id lookup
#   scDesign3/SymSim: fit_key lookup, reconstructed via the EXACT SAME apply()
#     call used during ground-truth extraction (column set + order matters --
#     this reproduces the padded string format already in each manifest)
#
# Rows with no matching ground truth (e.g. Splatter's 6 documented
# exclusions) are SKIPPED, not written -- same principle as the ground-
# truth extraction stage: never produce a plausible-looking artifact for
# data that hasn't been verified.
#
# Checkpointed, periodic progress logging, per-row error handling,
# persistent warning capture, and a live disk-space safety halt (every 200
# rows) -- all built in from the start this time.
# =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

pg <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

fit_cols_scd <- c("depth","separability","n_cells","batch")
pg$fit_key_scd <- apply(pg[,fit_cols_scd], 1, function(r) paste(r, collapse="_"))
fit_cols_sym <- c("separability","n_cells","batch","dropout","depth")
pg$fit_key_sym <- apply(pg[,fit_cols_sym], 1, function(r) paste(r, collapse="_"))

scd_manifest <- read.csv("data/simulated/ground_truth/scdesign3_manifest.csv", stringsAsFactors=FALSE)
sym_manifest <- read.csv("data/simulated/ground_truth/symsim_manifest.csv", stringsAsFactors=FALSE)

OUT_ROOT <- "data/simulated/sce"
LOG_DIR <- "logs"
MIN_FREE_GB <- 5
dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)

subset_env <- Sys.getenv("SCE_RUN_ID_SUBSET", unset="")
sim_env    <- Sys.getenv("SCE_SIMULATORS", unset="splatter,scdesign3,symsim")
sims <- strsplit(sim_env, ",")[[1]]

check_free_gb <- function(path=OUT_ROOT) {
  out <- tryCatch(system(sprintf("df -BG --output=avail %s | tail -1", path), intern=TRUE),
                   error=function(e) NA)
  suppressWarnings(as.numeric(gsub("[^0-9]", "", out)))
}

convert_one <- function(sim, run_id) {
  out_dir <- file.path(OUT_ROOT, sim)
  out_path <- file.path(out_dir, sprintf("%s_sce_run_%05d.rds", sim, run_id))
  if (file.exists(out_path)) return(list(status="skip", path=out_path))

  prod_path <- sprintf("data/simulated/%s/%s_run_%05d.rds", sim, sim, run_id)
  if (!file.exists(prod_path)) return(list(status="missing_production", path=NA))

  prod <- readRDS(prod_path)
  counts <- prod$counts
  if (sim == "symsim" && is.null(rownames(counts))) {
    rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
  }

  gt <- NULL
  gt_source <- NA
  if (sim == "splatter") {
    gt_path <- sprintf("data/simulated/ground_truth/splatter/splatter_truth_run_%05d.rds", run_id)
    if (file.exists(gt_path)) gt <- readRDS(gt_path)
  } else if (sim == "scdesign3") {
    fk <- pg$fit_key_scd[pg$run_id==run_id]
    match_path <- scd_manifest$output_path[scd_manifest$fit_key==fk]
    if (length(match_path)==1 && file.exists(match_path)) gt <- readRDS(match_path)
  } else if (sim == "symsim") {
    fk <- pg$fit_key_sym[pg$run_id==run_id]
    match_path <- sym_manifest$output_path[sym_manifest$fit_key==fk]
    if (length(match_path)==1 && file.exists(match_path)) gt <- readRDS(match_path)
  }
  if (is.null(gt)) return(list(status="no_ground_truth", path=NA))

  sce <- SingleCellExperiment(assays=list(counts=counts), colData=DataFrame(prod$cell_meta))
  metadata(sce)$run_params <- prod$run_params
  metadata(sce)$true_group_means <- gt$true_group_means
  metadata(sce)$ground_truth_source <- if (!is.null(gt$source)) gt$source else gt$method_note
  metadata(sce)$simulator <- sim

  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  tryCatch({
    saveRDS(sce, out_path, compress=TRUE)
  }, error = function(e) {
    if (file.exists(out_path)) file.remove(out_path)  # never leave a partial file behind
    stop(e)
  })
  list(status="done", path=out_path, dim=paste(dim(sce), collapse="x"),
       size_mb=file.info(out_path)$size/1e6)
}

manifest_path <- file.path(OUT_ROOT, "sce_conversion_manifest.csv")
dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(manifest_path)) {
  write.csv(data.frame(matrix(nrow=0,ncol=4,dimnames=list(NULL,c("simulator","run_id","status","output_path")))),
            manifest_path, row.names=FALSE)
}
warn_con <- file(file.path(LOG_DIR, "convert_to_sce_warnings.log"), open="a")

for (sim in sims) {
  ids <- pg$run_id
  if (nzchar(subset_env)) ids <- as.integer(strsplit(subset_env, ",")[[1]])

  cat(sprintf("\n=== %s: %d candidate rows ===\n", sim, length(ids)))
  t0 <- Sys.time(); n_done <- 0; n_skip <- 0; n_nogt <- 0; n_err <- 0

  for (i in seq_along(ids)) {
    run_id <- ids[i]

    if (i %% 200 == 0) {
      free_gb <- check_free_gb()
      if (!is.na(free_gb) && free_gb < MIN_FREE_GB) {
        msg <- sprintf("HALT: free space %sGB < %sGB threshold, stopping safely", free_gb, MIN_FREE_GB)
        cat(msg, "\n"); writeLines(msg, warn_con); flush(warn_con)
        stop(msg)
      }
    }

    res <- tryCatch(convert_one(sim, run_id), error=function(e) {
      msg <- sprintf("[%s] %s run_id=%d ERROR: %s", format(Sys.time()), sim, run_id, conditionMessage(e))
      writeLines(msg, warn_con); flush(warn_con)
      list(status="error", path=NA)
    })

    if (res$status %in% c("done","skip")) n_done <- n_done + (res$status=="done")
    if (res$status=="skip") n_skip <- n_skip + 1
    if (res$status=="no_ground_truth") n_nogt <- n_nogt + 1
    if (res$status=="error") n_err <- n_err + 1

    if (res$status != "skip") {
      write.table(data.frame(simulator=sim, run_id=run_id, status=res$status, output_path=res$path),
                  manifest_path, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }

    if (i %% 500 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      cat(sprintf("[%s] %s: %d/%d  done=%d skip=%d no_gt=%d err=%d  elapsed=%.1fmin\n",
                  format(Sys.time()), sim, i, length(ids), n_done, n_skip, n_nogt, n_err, elapsed/60))
    }
  }
  cat(sprintf("%s complete: done=%d skip=%d no_gt=%d err=%d\n", sim, n_done, n_skip, n_nogt, n_err))
}
close(warn_con)
