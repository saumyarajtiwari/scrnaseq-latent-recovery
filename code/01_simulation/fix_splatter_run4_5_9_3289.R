# =============================================================================
# fix_splatter_run4_5_9_3289.R
#
# Regenerates 4 Splatter main-grid files (run_id 4, 5, 9, 3289) confirmed
# corrupted: all four were the entire output of the FIRST invocation of
# simulate_splatter.R (2026-06-23 11:39:07-21), under BiocParallel's
# MulticoreParam, before it was replaced with parallel::mclapply (commit
# eaa5bf87, same day 16:10:03) specifically because of a parallel-backend
# defect. BiocParallel manages its own per-worker RNG substreams, which
# can override/interact with a seed argument passed into the dispatched
# function -- plausibly producing draws inconsistent with what the same
# seeded call produces under mclapply or a single-threaded call.
#
# Discovered during Step 1.8 ground-truth extraction validation (a 12-row
# test batch flagged run_id=4,5 with group_match=FALSE; 82-85% of cell-
# group assignments differed from a fresh regen -- the signature of an
# independent random draw, not a subtle parameter difference).
#
# Investigation ruled out the initially-suspected sparsity-key-lookup fix
# (5e59215b): same-sparsity-label files generated minutes later, under the
# corrected mclapply engine but still BEFORE that lookup fix landed, are
# confirmed correct -- isolating the parallel-backend switch as the actual
# causal event, not the lookup fix.
#
# Scope confirmed via direct regen-and-compare against all 128 files
# sharing the same early mtime window: only these 4 (0.037% of the
# 10,935-row main grid) are affected -- exactly the files produced before
# the 271-minute gap ending in the parallel-backend fix.
#
# param_grid.csv values for these 4 rows were always correct; only the
# resulting .rds content was wrong. No param_grid.csv change needed.
#
# Old (corrupted) files backed up, not deleted, before overwrite. Backup
# dir is local-only: data/simulated/**/*.rds is gitignored, so neither
# the original nor fixed .rds files were ever git-tracked -- only this
# fix script is.
# =============================================================================

suppressPackageStartupMessages({
  library(splatter)
  library(SingleCellExperiment)
  library(Matrix)
})

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

TARGET_IDS <- c(4, 5, 9, 3289)
N_GENES <- 10000L
OUT_DIR <- "data/simulated/splatter"

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
BACKUP_DIR <- sprintf("data/simulated/splatter_run4_5_9_3289_v1_INVALID_backup_%s", TIMESTAMP)
dir.create(BACKUP_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== Splatter corrupted-file fix: run_id", paste(TARGET_IDS, collapse=", "), "===\n")
cat("Backup dir:", BACKUP_DIR, "\n\n")

build_args <- function(row) {
  sparsity_p     <- splatter_sparsity[[sprintf("%.2f", as.numeric(row$sparsity))]]
  depth_p        <- splatter_depth[[as.character(row$depth)]]
  dropout_p      <- splatter_dropout[[row$dropout]]
  separability_p <- splatter_separability[[row$separability]]
  n_cells        <- as.integer(row$n_cells)
  is_null        <- isTRUE(as.logical(row$is_null_control)) || as.integer(row$n_groups) == 1L

  if (row$batch == "none") {
    batch_cells <- n_cells; bfl <- 0.1; bfs <- 0.1
  } else if (row$batch == "simple") {
    h <- floor(n_cells / 2L)
    batch_cells <- c(h, n_cells - h); bfl <- 0.1; bfs <- 0.1
  } else {
    b1 <- round(n_cells * 0.50); b2 <- round(n_cells * 0.30)
    batch_cells <- c(b1, b2, n_cells - b1 - b2); bfl <- 0.2; bfs <- 0.2
  }

  args <- list(nGenes=N_GENES, batchCells=batch_cells, bcv.common=sparsity_p$bcv.common,
               lib.loc=depth_p$lib.loc, lib.scale=depth_p$lib.scale,
               dropout.type=dropout_p$dropout.type, batch.facLoc=bfl, batch.facScale=bfs,
               verbose=FALSE, seed=as.integer(row$run_id))
  if (!is.null(dropout_p$dropout.mid)) args$dropout.mid <- dropout_p$dropout.mid

  if (is_null) {
    args$method <- "single"
  } else {
    args$method <- "groups"; args$group.prob <- rep(1/5, 5)
    args$de.prob <- separability_p$de.prob
    args$de.facLoc <- separability_p$de.facLoc
    args$de.facScale <- separability_p$de.facScale
  }
  list(args=args, is_null=is_null, sparsity_p=sparsity_p, depth_p=depth_p, n_cells=n_cells)
}

regenerate_one <- function(run_id) {
  cat(sprintf("--- run_id=%05d ---\n", run_id))
  row <- param_grid[param_grid$run_id == run_id, ]
  stopifnot(nrow(row) == 1)

  out_path <- file.path(OUT_DIR, sprintf("splatter_run_%05d.rds", run_id))
  backup_path <- file.path(BACKUP_DIR, sprintf("splatter_run_%05d.rds", run_id))
  stopifnot(file.exists(out_path))
  file.copy(out_path, backup_path, overwrite = FALSE)
  cat("Backed up old file to:", backup_path, "\n")

  b <- build_args(row)
  sim  <- do.call(splatSimulate, b$args)
  cnts <- as(counts(sim), "dgCMatrix")
  actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)

  true_group <- if ("Group" %in% names(colData(sim))) as.character(sim$Group) else rep("Group1", ncol(sim))
  batch_id <- if ("Batch" %in% names(colData(sim))) {
    as.integer(sub("Batch", "", as.character(sim$Batch)))
  } else rep(1L, ncol(sim))

  cell_meta <- data.frame(run_id=run_id, cell_id=colnames(sim), true_group=true_group,
                           batch_id=batch_id, achieved_sparsity=actual_sparsity, stringsAsFactors=FALSE)

  run_params <- list(run_id=run_id, sparsity_label=as.character(row$sparsity),
                      depth_label=as.character(row$depth), dropout=row$dropout,
                      separability=row$separability, n_cells=b$n_cells,
                      n_groups=if (b$is_null) 1L else 5L, batch=row$batch,
                      gene_strategy=row$gene_strategy, clipping=row$clipping,
                      is_null_control=b$is_null, bcv_common=b$sparsity_p$bcv.common,
                      lib_loc=b$depth_p$lib.loc, actual_sparsity=actual_sparsity)

  saveRDS(list(counts=cnts, cell_meta=cell_meta, run_params=run_params), file=out_path, compress=TRUE)
  cat(sprintf("Regenerated. n_groups=%d  achieved_sparsity=%.4f\n", length(unique(true_group)), actual_sparsity))
  rm(sim); gc(verbose=FALSE)
}

for (rid in TARGET_IDS) regenerate_one(rid)

cat("\n=== Independent re-verification: THIRD fresh call vs newly-written file ===\n")
verify_one <- function(run_id) {
  row <- param_grid[param_grid$run_id == run_id, ]
  b <- build_args(row)
  sim <- do.call(splatSimulate, b$args)
  fresh_group  <- if ("Group" %in% names(colData(sim))) as.character(sim$Group) else rep("Group1", ncol(sim))
  fresh_counts <- as(counts(sim), "dgCMatrix")

  new_file <- readRDS(sprintf("data/simulated/splatter/splatter_run_%05d.rds", run_id))
  group_ok  <- identical(fresh_group, new_file$cell_meta$true_group)
  counts_ok <- identical(fresh_counts, new_file$counts)
  rm(sim); gc(verbose=FALSE)
  cat(sprintf("run_id=%05d  group_identical=%s  counts_identical=%s\n", run_id, group_ok, counts_ok))
  c(group_ok, counts_ok)
}
all_ok <- sapply(TARGET_IDS, verify_one)
cat(sprintf("\nALL VERIFIED: %s\n", all(all_ok)))
