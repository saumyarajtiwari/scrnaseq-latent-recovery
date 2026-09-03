# =============================================================================
# extract_ground_truth_splatter.R
#
# Extracts per-run_id "true" biological signal for all 10,935 main-grid
# Splatter runs: group-mean expression from the TrueCounts assay (confirmed
# present under both method="groups" and method="single"; confirmed to
# correctly carry DE signal via direct DEFacGroup correlation check).
# No fit-key batching exists for Splatter (seed=run_id unique per row).
#
# Validated per row against the already-verified production file:
#   - group_match: exact identity of true_group label vector
#   - value_match:
#       dropout=="none"       -> TrueCounts and production counts must be
#                                 EXACTLY identical (confirmed directly)
#       dropout in {low,high} -> production counts must equal TrueCounts
#                                 exactly at every position production kept
#                                 nonzero (dropout confirmed to be pure
#                                 zero-masking, never modifies surviving
#                                 values), plus nnz/sum bounds
#
# ROWS WITH A CONFIRMED MISMATCH ARE NOT WRITTEN. A plausible-looking
# ground-truth file paired with a different production count matrix is a
# real correctness risk if anything downstream doesn't specifically check
# the match flags -- so such rows are logged to splatter_unresolved.csv
# instead, with no .rds written. Checkpointing (skip-if-exists) means a
# later re-run automatically picks up and correctly processes any such
# row once its underlying production file is fixed -- no special handling
# needed. (This design was itself directly motivated by a real finding:
# run_id 4/5/9/3289 were confirmed corrupted -- the entire output of an
# abandoned first invocation of simulate_splatter.R under a since-replaced
# parallel backend -- and have since been regenerated and fixed.)
#
# Checkpointed (skip-if-exists). Progress logged every 50 attempted rows
# (written+unresolved+errors) with timestamps. Warnings captured
# persistently, one log file per worker. Per-row errors are caught,
# logged, and skipped -- one anomalous row cannot halt an unattended
# multi-hour run.
#
# sim (the full SCE, ~1.6GB resting / ~3GB transient at worst-case scale,
# 4 unused dense assays) is freed immediately after TrueCounts/Group are
# extracted, not held through validation+save.
#
# Optional env vars:
#   GT_RUN_ID_SUBSET=1,31,16     restrict to specific run_ids (testing)
#   GT_WORKER_COUNT=2            total parallel workers for a real split
#   GT_WORKER_INDEX=0            this worker's 0-indexed slot
#   GT_PROD_DIR_OVERRIDE=<path>  redirect the "production" comparison dir
#                                (default: data/simulated/splatter) --
#                                used for testing the skip-and-log path
#                                against a known-bad backup copy without
#                                touching real production data
#
# Read-only w.r.t. param_grid.csv, param_dict.R, and the production
# directory (real or overridden). Writes only to:
#   data/simulated/ground_truth/splatter/splatter_truth_run_NNNNN.rds
#   data/simulated/ground_truth/splatter_manifest.csv
#   data/simulated/ground_truth/splatter_unresolved.csv
#   logs/extract_ground_truth_splatter_warnings[_wN].log
# =============================================================================

suppressPackageStartupMessages({
  library(splatter)
  library(SingleCellExperiment)
  library(Matrix)
})

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)
main_grid  <- param_grid[param_grid$is_null_control == FALSE, ]

subset_env <- Sys.getenv("GT_RUN_ID_SUBSET", unset = "")
if (nzchar(subset_env)) {
  subset_ids <- as.integer(strsplit(subset_env, ",")[[1]])
  main_grid <- main_grid[main_grid$run_id %in% subset_ids, ]
}

worker_count <- as.integer(Sys.getenv("GT_WORKER_COUNT", unset = "1"))
worker_index <- as.integer(Sys.getenv("GT_WORKER_INDEX", unset = "0"))
if (worker_count > 1) {
  keep <- (seq_len(nrow(main_grid)) %% worker_count) == (worker_index %% worker_count)
  main_grid <- main_grid[keep, ]
}

N_GENES  <- 10000L
OUT_DIR  <- "data/simulated/ground_truth/splatter"
PROD_DIR <- Sys.getenv("GT_PROD_DIR_OVERRIDE", unset = "data/simulated/splatter")
LOG_DIR  <- "logs"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)

worker_tag <- if (worker_count > 1) sprintf("_w%d", worker_index) else ""
warn_log <- file.path(LOG_DIR, sprintf("extract_ground_truth_splatter_warnings%s.log", worker_tag))
warn_con <- file(warn_log, open="a")

cat("=== Splatter Ground-Truth Extraction ===\n")
if (nzchar(subset_env)) cat(sprintf("*** GT_RUN_ID_SUBSET active: %s ***\n", subset_env))
if (worker_count > 1) cat(sprintf("*** Worker %d of %d (0-indexed) ***\n", worker_index, worker_count))
if (PROD_DIR != "data/simulated/splatter") cat(sprintf("*** GT_PROD_DIR_OVERRIDE active: %s ***\n", PROD_DIR))
cat(sprintf("Main-grid rows   : %d\n", nrow(main_grid)))
cat(sprintf("Production dir   : %s\n", PROD_DIR))
cat(sprintf("Started          : %s\n", format(Sys.time())))
cat(sprintf("Warning log      : %s\n\n", warn_log))

manifest_path <- "data/simulated/ground_truth/splatter_manifest.csv"
manifest_cols <- c("run_id","n_groups","group_match","value_match","output_path")
if (!file.exists(manifest_path)) {
  write.csv(data.frame(matrix(nrow=0, ncol=length(manifest_cols),
                               dimnames=list(NULL, manifest_cols))),
            manifest_path, row.names=FALSE)
}

unresolved_path <- "data/simulated/ground_truth/splatter_unresolved.csv"
unresolved_cols <- c("run_id","group_match","value_match","timestamp","note")
if (!file.exists(unresolved_path)) {
  write.csv(data.frame(matrix(nrow=0, ncol=length(unresolved_cols),
                               dimnames=list(NULL, unresolved_cols))),
            unresolved_path, row.names=FALSE)
}

t_start <- Sys.time()
n_mismatch <- 0L
n_error <- 0L
n_done_this_run <- 0L

for (i in seq_len(nrow(main_grid))) {
  row    <- main_grid[i, ]
  run_id <- as.integer(row$run_id)
  out_path <- file.path(OUT_DIR, sprintf("splatter_truth_run_%05d.rds", run_id))

  if (file.exists(out_path)) next

  result <- tryCatch({

    sparsity_p     <- splatter_sparsity[[sprintf("%.2f", as.numeric(row$sparsity))]]
    depth_p        <- splatter_depth[[as.character(row$depth)]]
    dropout_p      <- splatter_dropout[[row$dropout]]
    separability_p <- splatter_separability[[row$separability]]
    n_cells        <- as.integer(row$n_cells)

    if (row$batch == "none") {
      batch_cells <- n_cells; batch_fac_loc <- 0.1; batch_fac_scale <- 0.1
    } else if (row$batch == "simple") {
      h <- floor(n_cells / 2L)
      batch_cells <- c(h, n_cells - h); batch_fac_loc <- 0.1; batch_fac_scale <- 0.1
    } else {
      b1 <- round(n_cells * 0.50); b2 <- round(n_cells * 0.30)
      batch_cells <- c(b1, b2, n_cells - b1 - b2); batch_fac_loc <- 0.2; batch_fac_scale <- 0.2
    }

    args <- list(
      nGenes = N_GENES, batchCells = batch_cells,
      bcv.common = sparsity_p$bcv.common,
      lib.loc = depth_p$lib.loc, lib.scale = depth_p$lib.scale,
      dropout.type = dropout_p$dropout.type,
      batch.facLoc = batch_fac_loc, batch.facScale = batch_fac_scale,
      verbose = FALSE, seed = run_id,
      method = "groups", group.prob = rep(1/5, 5),
      de.prob = separability_p$de.prob,
      de.facLoc = separability_p$de.facLoc,
      de.facScale = separability_p$de.facScale
    )
    if (!is.null(dropout_p$dropout.mid)) args$dropout.mid <- dropout_p$dropout.mid

    withCallingHandlers({

      sim <- do.call(splatSimulate, args)
      tc  <- as(assay(sim, "TrueCounts"), "dgCMatrix")
      true_group <- if ("Group" %in% names(colData(sim))) as.character(sim$Group) else rep("Group1", ncol(sim))
      rm(sim); gc(verbose=FALSE, full=TRUE)

      rownames(tc) <- paste0("Gene", seq_len(nrow(tc)))
      group_levels <- sort(unique(true_group))
      true_group_means <- vapply(group_levels, function(g)
        Matrix::rowMeans(as.matrix(tc[, true_group==g, drop=FALSE])),
        FUN.VALUE = numeric(nrow(tc)))
      colnames(true_group_means) <- group_levels

      prod_path <- file.path(PROD_DIR, sprintf("splatter_run_%05d.rds", run_id))
      group_match <- value_match <- NA

      if (file.exists(prod_path)) {
        prod <- readRDS(prod_path)
        group_match <- identical(true_group, prod$cell_meta$true_group)

        if (row$dropout == "none") {
          diff <- tc - prod$counts
          value_match <- (Matrix::nnzero(diff) == 0)
        } else {
          prod_t <- as(prod$counts, "TsparseMatrix")
          tc_vals <- tc[cbind(prod_t@i + 1L, prod_t@j + 1L)]
          exact_at_kept <- all(tc_vals == prod_t@x)
          nnz_ok <- Matrix::nnzero(tc) >= Matrix::nnzero(prod$counts)
          sum_ok <- sum(tc) >= sum(prod$counts)
          value_match <- exact_at_kept && nnz_ok && sum_ok
        }
      }

      verified <- isTRUE(group_match) && isTRUE(value_match)

      if (verified) {
        saveRDS(list(
          true_group_means = true_group_means,
          run_id = run_id,
          group_levels = group_levels,
          source = "splatter_TrueCounts_assay",
          group_match = group_match,
          value_match = value_match,
          method_note = "Group means of TrueCounts assay (pre-dropout). Validated against already-verified production: exact identity when dropout==none, exact match at kept-nonzero positions plus nnz/sum bounds otherwise."
        ), out_path, compress=TRUE)

        manifest_row <- data.frame(run_id=run_id, n_groups=length(group_levels),
                                    group_match=group_match, value_match=value_match,
                                    output_path=out_path)
        write.table(manifest_row, manifest_path, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
      } else {
        n_mismatch <<- n_mismatch + 1L
        msg <- sprintf("UNRESOLVED run_id=%05d: group_match=%s value_match=%s -- NOT written, logged only",
                        run_id, group_match, value_match)
        cat(msg, "\n"); writeLines(msg, warn_con); flush(warn_con)

        unresolved_row <- data.frame(run_id=run_id, group_match=group_match, value_match=value_match,
                                      timestamp=format(Sys.time()),
                                      note="not written -- production mismatch, pending production fix")
        write.table(unresolved_row, unresolved_path, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
      }

      rm(tc)
      verified

    }, warning = function(w) {
      msg <- sprintf("[%s] run_id=%05d WARNING: %s", format(Sys.time()), run_id, conditionMessage(w))
      writeLines(msg, warn_con); flush(warn_con)
      invokeRestart("muffleWarning")
    })

  }, error = function(e) {
    n_error <<- n_error + 1L
    msg <- sprintf("[%s] run_id=%05d ERROR (row skipped, no output written): %s",
                    format(Sys.time()), run_id, conditionMessage(e))
    cat(msg, "\n"); writeLines(msg, warn_con); flush(warn_con)
    FALSE
  })

  if (isTRUE(result)) n_done_this_run <- n_done_this_run + 1L

  if ((n_done_this_run + n_error + n_mismatch) %% 50 == 0 && (n_done_this_run + n_error + n_mismatch) > 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
    cat(sprintf("[%s] %d written + %d unresolved + %d errors (i=%d/%d)  elapsed=%.1fmin  rate=%.2fs/row\n",
                format(Sys.time()), n_done_this_run, n_mismatch, n_error, i, nrow(main_grid),
                elapsed/60, elapsed/(n_done_this_run+n_error+n_mismatch)))
  }
}

close(warn_con)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
cat(sprintf("\n=== DONE: %d written, %d unresolved (not written), %d errors, this run in %.1f min ===\n",
            n_done_this_run, n_mismatch, n_error, elapsed/60))
