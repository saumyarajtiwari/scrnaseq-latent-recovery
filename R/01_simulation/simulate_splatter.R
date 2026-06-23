# =============================================================================
# simulate_splatter.R
#
# Runs all 10,940 Splatter simulations from param_grid.csv.
# One .rds output per run: sparse count matrix + cell metadata + run parameters.
#
# Key decisions:
#   nGenes       = 10,000 (fixed throughout)
#   N_WORKERS    = 10     (parallel::mclapply, fork-based)
#   Checkpointing: re-running skips already-completed files
#   gene_strategy and clipping: preprocessing-stage instructions only —
#     recorded in run_params but do not affect the simulation
#
# Output: data/simulated/splatter/splatter_run_NNNNN.rds
# Log:    logs/splatter_run.log
# =============================================================================

suppressPackageStartupMessages({
  library(splatter)
  library(Matrix)
  library(parallel)
})

# -----------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors = FALSE)

N_GENES   <- 10000L
N_WORKERS <- 10L
OUT_DIR   <- "data/simulated/splatter"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Splatter Simulation ===\n")
cat(sprintf("Total runs       : %d\n",   nrow(param_grid)))
cat(sprintf("Genes per run    : %d\n",   N_GENES))
cat(sprintf("Parallel workers : %d\n",   N_WORKERS))
cat(sprintf("Output directory : %s\n",   OUT_DIR))
cat(sprintf("Started          : %s\n\n", format(Sys.time())))

# -----------------------------------------------------------------------
# PER-RUN FUNCTION
# -----------------------------------------------------------------------

run_one_splatter <- function(i) {

  row      <- param_grid[i, ]
  run_id   <- row$run_id
  out_path <- file.path(OUT_DIR, sprintf("splatter_run_%05d.rds", run_id))

  # Checkpointing: skip completed runs
  if (file.exists(out_path)) {
    cat(sprintf("[%05d] SKIP\n", run_id))
    return(run_id)
  }

  # -- Parameters from param_dict --
  sparsity_p     <- splatter_sparsity[[sprintf("%.2f", as.numeric(row$sparsity))]]
  depth_p        <- splatter_depth[[as.character(row$depth)]]
  dropout_p      <- splatter_dropout[[row$dropout]]
  separability_p <- splatter_separability[[row$separability]]
  n_cells        <- as.integer(row$n_cells)

  # -- Batch configuration --
  if (row$batch == "none") {
    batch_cells     <- n_cells
    batch_fac_loc   <- 0.1
    batch_fac_scale <- 0.1
  } else if (row$batch == "simple") {
    h <- floor(n_cells / 2L)
    batch_cells     <- c(h, n_cells - h)          # 50/50 split
    batch_fac_loc   <- 0.1
    batch_fac_scale <- 0.1
  } else {                                          # complex: 50/30/20 split
    b1 <- round(n_cells * 0.50)
    b2 <- round(n_cells * 0.30)
    batch_cells     <- c(b1, b2, n_cells - b1 - b2)
    batch_fac_loc   <- 0.2
    batch_fac_scale <- 0.2
  }

  # -- Null control vs multi-group --
  is_null <- isTRUE(as.logical(row$is_null_control)) || as.integer(row$n_groups) == 1L

  # -- Assemble splatSimulate arguments --
  args <- list(
    nGenes         = N_GENES,
    batchCells     = batch_cells,
    bcv.common     = sparsity_p$bcv.common,
    lib.loc        = depth_p$lib.loc,
    lib.scale      = depth_p$lib.scale,
    dropout.type   = dropout_p$dropout.type,
    batch.facLoc   = batch_fac_loc,
    batch.facScale = batch_fac_scale,
    verbose        = FALSE,
    seed           = run_id
  )

  if (!is.null(dropout_p$dropout.mid)) {
    args$dropout.mid <- dropout_p$dropout.mid
  }

  if (is_null) {
    args$method <- "single"
  } else {
    args$method      <- "groups"
    args$group.prob  <- rep(1/5, 5)
    args$de.prob     <- separability_p$de.prob
    args$de.facLoc   <- separability_p$de.facLoc
    args$de.facScale <- separability_p$de.facScale
  }

  # -- Simulate and save --
  tryCatch({

    sim  <- do.call(splatSimulate, args)
    cnts <- as(counts(sim), "dgCMatrix")

    actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)

    true_group <- if ("Group" %in% names(colData(sim))) {
      as.character(sim$Group)
    } else {
      rep("Group1", ncol(sim))
    }

    batch_id <- if ("Batch" %in% names(colData(sim))) {
      as.integer(sub("Batch", "", as.character(sim$Batch)))
    } else {
      rep(1L, ncol(sim))
    }

    cell_meta <- data.frame(
      run_id            = run_id,
      cell_id           = colnames(sim),
      true_group        = true_group,
      batch_id          = batch_id,
      achieved_sparsity = actual_sparsity,
      stringsAsFactors  = FALSE
    )

    run_params <- list(
      run_id          = run_id,
      sparsity_label  = as.character(row$sparsity),
      depth_label     = as.character(row$depth),
      dropout         = row$dropout,
      separability    = row$separability,
      n_cells         = n_cells,
      n_groups        = if (is_null) 1L else 5L,
      batch           = row$batch,
      gene_strategy   = row$gene_strategy,   # preprocessing instruction only
      clipping        = row$clipping,         # preprocessing instruction only
      is_null_control = is_null,
      bcv_common      = sparsity_p$bcv.common,
      lib_loc         = depth_p$lib.loc,
      actual_sparsity = actual_sparsity
    )

    saveRDS(
      object   = list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
      file     = out_path,
      compress = TRUE
    )

    cat(sprintf("[%05d] DONE   sparsity=%.4f  cells=%d  %s\n",
                run_id, actual_sparsity, ncol(sim),
                format(Sys.time(), "%H:%M:%S")))

    return(run_id)

  }, error = function(e) {
    cat(sprintf("[%05d] ERROR  %s\n", run_id, conditionMessage(e)))
    return(-run_id)
  })
}

# -----------------------------------------------------------------------
# PARALLEL EXECUTION
# -----------------------------------------------------------------------

cat(sprintf("Launching %d workers via parallel::mclapply ...\n\n", N_WORKERS))

mclapply(
  X              = seq_len(nrow(param_grid)),
  FUN            = run_one_splatter,
  mc.cores       = N_WORKERS,
  mc.preschedule = FALSE
)

# -----------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------

n_done   <- length(list.files(OUT_DIR, pattern = "^splatter_run_.*\\.rds$"))
n_target <- nrow(param_grid)

cat(sprintf("\n=== Run complete: %s ===\n", format(Sys.time())))
cat(sprintf("Completed : %d / %d\n", n_done, n_target))
if (n_done < n_target) {
  cat(sprintf("Failed    : %d  — re-run this script to retry\n", n_target - n_done))
}
