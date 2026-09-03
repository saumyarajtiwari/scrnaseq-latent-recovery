# =============================================================================
# simulate_symsim.R
#
# Runs all 10,940 SymSim simulations from param_grid.csv.
# One .rds output per run: sparse count matrix + cell metadata + run parameters.
#
# Strategy: simulate-once-save-many
#   Unique SymSim calls: separability(4) x n_cells(3) x batch(3) x
#                        dropout(3) x depth(3) = 244 unique simulations
#   Sparsity label does not affect simulation (Sigma fixed at 0.4)
#   gene_strategy and clipping are preprocessing-stage labels only
#   Each unique call saves 45 files (main) or 5 files (null)
#
# Reproducibility:
#   SimulateTrueCounts  -> randseed = run_id of first row in fit group
#   True2ObservedCounts -> set.seed(run_id of first row) before each call
#
# Batch: independent sub-simulations per batch, concatenated
#   none:    seed = run_id
#   simple:  seeds = run_id+10000, run_id+20000
#   complex: seeds = run_id+10000, run_id+20000, run_id+30000
#
# Output: data/simulated/symsim/symsim_run_NNNNN.rds
# Log:    logs/symsim_run.log
# =============================================================================

suppressPackageStartupMessages({
  library(SymSim)
  library(Matrix)
  library(parallel)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

N_GENES   <- 2000L
N_WORKERS <- 2L
OUT_DIR   <- "data/simulated/symsim"

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== SymSim Simulation ===\n")
cat(sprintf("Total runs       : %d\n",   nrow(param_grid)))
cat(sprintf("Genes per run    : %d\n",   N_GENES))
cat(sprintf("Parallel workers : %d\n",   N_WORKERS))
cat(sprintf("Output directory : %s\n",   OUT_DIR))
cat(sprintf("Started          : %s\n\n", format(Sys.time())))

# =============================================================================
# IDENTIFY UNIQUE FIT GROUPS
# Sparsity, gene_strategy, clipping do not affect SymSim output
# =============================================================================

fit_cols <- c("separability", "n_cells", "batch", "dropout", "depth")
param_grid$fit_key <- apply(param_grid[, fit_cols], 1,
                            function(r) paste(r, collapse="_"))
unique_fits <- unique(param_grid$fit_key)
cat(sprintf("Unique SymSim calls: %d\n\n", length(unique_fits)))

# =============================================================================
# PER-FIT FUNCTION
# =============================================================================

run_one_fit <- function(fk) {

  rows   <- param_grid[param_grid$fit_key == fk, ]
  row1   <- rows[1, ]
  run_id <- as.integer(row1$run_id)

  # Checkpointing — skip if all output files for this fit exist
  out_paths <- file.path(OUT_DIR,
    sprintf("symsim_run_%05d.rds", rows$run_id))
  if (all(file.exists(out_paths))) {
    cat(sprintf("[fit %s] SKIP\n", fk))
    return(invisible(NULL))
  }

  tryCatch({

    # ----------------------------------------------------------------
    # Parameters
    # ----------------------------------------------------------------
    sep_p     <- symsim_separability[[row1$separability]]
    dropout_p <- symsim_dropout[[row1$dropout]]
    depth_p   <- symsim_depth[[row1$dropout]][[as.character(row1$depth)]]
    n_cells   <- as.integer(row1$n_cells)

    is_null  <- isTRUE(as.logical(row1$is_null_control)) ||
                row1$separability == "null"
    evf_type <- sep_p$evf_type
    n_de_evf <- if (is_null) 8L else as.integer(sep_p$n_de_evf)

    # ----------------------------------------------------------------
    # Batch configuration
    # ----------------------------------------------------------------
    if (row1$batch == "none") {
      batch_sizes <- n_cells
      batch_seeds <- run_id
    } else if (row1$batch == "simple") {
      b1 <- floor(n_cells / 2L)
      batch_sizes <- c(b1, n_cells - b1)
      batch_seeds <- c(run_id + 10000L, run_id + 20000L)
    } else {
      b1 <- round(n_cells * 0.50)
      b2 <- round(n_cells * 0.30)
      batch_sizes <- c(b1, b2, n_cells - b1 - b2)
      batch_seeds <- c(run_id + 10000L, run_id + 20000L, run_id + 30000L)
    }

    # ----------------------------------------------------------------
    # Simulate each batch and concatenate
    # ----------------------------------------------------------------
    all_counts <- vector("list", length(batch_sizes))
    all_pops   <- vector("list", length(batch_sizes))
    all_batch  <- vector("list", length(batch_sizes))

    for (b in seq_along(batch_sizes)) {
      bs    <- as.integer(batch_sizes[b])
      bseed <- as.integer(batch_seeds[b])
      min_pop <- if (is_null) bs else max(5L, floor(bs / 5L))

      tr <- SimulateTrueCounts(
        ncells_total = bs,
        min_popsize  = min_pop,
        i_minpop     = 1L,
        ngenes       = N_GENES,
        evf_type     = evf_type,
        phyla        = if (evf_type == "discrete") Phyla5() else NULL,
        randseed     = bseed,
        Sigma        = SYMSIM_SIGMA,
        n_de_evf     = n_de_evf,
        bimod        = 0,
        vary         = "s"
      )

      true_counts <- tr[["counts"]]
      cell_meta   <- tr[["cell_meta"]]

      set.seed(bseed)
      obs <- True2ObservedCounts(
        true_counts = true_counts,
        meta_cell   = cell_meta,
        protocol    = "UMI",
        alpha_mean  = dropout_p$alpha_mean,
        alpha_sd    = dropout_p$alpha_sd,
        gene_len    = rep(1000L, N_GENES),
        depth_mean  = depth_p$depth_mean,
        depth_sd    = depth_p$depth_sd
      )

      all_counts[[b]] <- obs[[1]]
      all_pops[[b]]   <- as.character(cell_meta$pop)
      all_batch[[b]]  <- rep(b, bs)
    }

    # Combine batches
    cnts            <- as(do.call(cbind, all_counts), "dgCMatrix")
    actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)
    n_total         <- ncol(cnts)
    true_group      <- unlist(all_pops)
    batch_id        <- as.integer(unlist(all_batch))

    # ----------------------------------------------------------------
    # Save one file per grid row sharing this fit
    # ----------------------------------------------------------------
    for (j in seq_len(nrow(rows))) {
      row_j    <- rows[j, ]
      out_path <- file.path(OUT_DIR,
        sprintf("symsim_run_%05d.rds", as.integer(row_j$run_id)))

      if (file.exists(out_path)) next

      cell_meta_out <- data.frame(
        run_id            = as.integer(row_j$run_id),
        cell_id           = paste0("cell_", seq_len(n_total)),
        true_group        = true_group,
        batch_id          = batch_id,
        achieved_sparsity = actual_sparsity,
        stringsAsFactors  = FALSE
      )

      run_params <- list(
        run_id          = as.integer(row_j$run_id),
        sparsity_label  = as.character(row_j$sparsity),
        depth_label     = as.character(row_j$depth),
        dropout         = row_j$dropout,
        separability    = row_j$separability,
        n_cells         = n_cells,
        n_cells_actual  = n_total,
        n_groups        = length(unique(true_group)),
        batch           = row_j$batch,
        gene_strategy   = row_j$gene_strategy,
        clipping        = row_j$clipping,
        is_null_control = is_null,
        sigma           = SYMSIM_SIGMA,
        n_de_evf        = n_de_evf,
        alpha_mean      = dropout_p$alpha_mean,
        depth_mean      = depth_p$depth_mean,
        actual_sparsity = actual_sparsity,
        fit_key         = fk
      )

      saveRDS(
        object   = list(counts=cnts, cell_meta=cell_meta_out,
                        run_params=run_params),
        file     = out_path,
        compress = TRUE
      )
    }

    cat(sprintf("[fit %s] DONE  sparsity=%.4f  cells=%d  groups=%d  files=%d  %s\n",
      fk, actual_sparsity, n_total, length(unique(true_group)),
      nrow(rows), format(Sys.time(), "%H:%M:%S")))

  }, error = function(e) {
    cat(sprintf("[fit %s] ERROR  %s\n", fk, conditionMessage(e)))
  })

  invisible(NULL)
}

# =============================================================================
# PARALLEL EXECUTION
# =============================================================================

cat(sprintf("Launching %d workers via parallel::mclapply ...\n\n", N_WORKERS))

mclapply(
  X              = unique_fits,
  FUN            = run_one_fit,
  mc.cores       = N_WORKERS,
  mc.preschedule = FALSE
)

# =============================================================================
# SUMMARY
# =============================================================================

n_done   <- length(list.files(OUT_DIR, pattern="^symsim_run_.*\\.rds$"))
n_target <- nrow(param_grid)

cat(sprintf("\n=== Run complete: %s ===\n", format(Sys.time())))
cat(sprintf("Completed : %d / %d\n", n_done, n_target))
if (n_done < n_target)
  cat(sprintf("Failed    : %d  — re-run to retry\n", n_target - n_done))
