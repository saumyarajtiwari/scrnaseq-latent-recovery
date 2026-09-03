# =============================================================================
# simulate_symsim.R  (REVISED — corrects a Step 1.7-discovered defect)
#
# CORRECTION NOTE:
#   The original version (preserved as simulate_symsim_v1_INVALID_sparsity_inert.R)
#   generated ONE count matrix per fit_key (separability x n_cells x batch x
#   dropout x depth — sparsity was never part of it, by explicit design: see
#   the original comment "Sparsity label does not affect simulation (Sigma
#   fixed at 0.4)"). That count matrix was saved unchanged across all 5
#   sparsity_label values sharing that fit_key. Result: sparsity_label had
#   ZERO effect on simulated data, despite being logged in run_params as if
#   it did. Discovered during Step 1.7 output validation (byte-identical
#   actual_sparsity across all 5 sparsity labels within every sampled group).
#
#   FIX: dropout (via alpha_mean/depth_mean, already genuinely working and
#   unchanged here) stays exactly as before. sparsity_label is now
#   implemented as explicit, calibrated post-hoc stochastic zero-masking
#   applied to the generated count matrix, using the same mechanism and
#   calibration methodology already validated for scDesign3 (see
#   simulate_scdesign3.R correction note):
#     - sparsity mask_p (0/.15/.35/.55/.75) from symsim_calib_sparsity.csv
#   Verified: strict monotonic sparsity ordering, zero degenerate all-zero
#   cells even at the worst-case combination tested (depth=500, dropout=high,
#   separability=high, smallest per-population size).
#
#   fit_key is unchanged (sparsity was never part of it) — same 244 unique
#   simulations as before, masking applied per sparsity_label when saving.
#
# Reproducibility:
#   SimulateTrueCounts  -> randseed = run_id of first row in fit group
#   True2ObservedCounts -> set.seed(run_id of first row) before each call
#   Post-hoc masking    -> set.seed(run_id * 1000 + sparsity_index), logged
#                          per file as mask_seed
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
sparsity_calib <- read.csv("data/simulated/symsim_calib_sparsity.csv", stringsAsFactors=FALSE)
mask_p_for_sparsity <- setNames(sparsity_calib$mask_p, as.character(sparsity_calib$sparsity_label))

N_GENES   <- 2000L
N_WORKERS <- 2L
OUT_DIR   <- "data/simulated/symsim"

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== SymSim Simulation (REVISED) ===\n")
cat(sprintf("Total runs       : %d\n",   nrow(param_grid)))
cat(sprintf("Genes per run    : %d\n",   N_GENES))
cat(sprintf("Parallel workers : %d\n",   N_WORKERS))
cat(sprintf("Output directory : %s\n",   OUT_DIR))
cat(sprintf("Started          : %s\n\n", format(Sys.time())))

# =============================================================================
# MASKING FUNCTION
# =============================================================================

apply_mask <- function(mat, p, seed) {
  if (p <= 0) return(mat)
  set.seed(seed)
  mat_t <- as(mat, "TsparseMatrix")
  n_nz <- length(mat_t@x)
  keep <- rbinom(n_nz, 1, 1 - p) == 1
  mat_t@x[!keep] <- 0
  drop0(as(mat_t, "CsparseMatrix"))
}

# =============================================================================
# IDENTIFY UNIQUE FIT GROUPS (unchanged — sparsity was never part of this)
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

  out_paths <- file.path(OUT_DIR,
    sprintf("symsim_run_%05d.rds", rows$run_id))
  if (all(file.exists(out_paths))) {
    cat(sprintf("[fit %s] SKIP\n", fk))
    return(invisible(NULL))
  }

  tryCatch({

    sep_p     <- symsim_separability[[row1$separability]]
    dropout_p <- symsim_dropout[[row1$dropout]]
    depth_p   <- symsim_depth[[row1$dropout]][[as.character(row1$depth)]]
    n_cells   <- as.integer(row1$n_cells)

    is_null  <- isTRUE(as.logical(row1$is_null_control)) ||
                row1$separability == "null"
    evf_type <- sep_p$evf_type
    n_de_evf <- if (is_null) 8L else as.integer(sep_p$n_de_evf)

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

    cnts_base       <- as(do.call(cbind, all_counts), "dgCMatrix")
    n_total         <- ncol(cnts_base)
    true_group      <- unlist(all_pops)
    batch_id        <- as.integer(unlist(all_batch))

    # ----------------------------------------------------------------
    # Apply sparsity masking once per distinct sparsity_label in this group
    # ----------------------------------------------------------------
    sparsity_labels_here <- unique(as.character(rows$sparsity))
    masked_cache <- list()
    for (k in seq_along(sparsity_labels_here)) {
      s_lbl <- sparsity_labels_here[k]
      p_sp  <- mask_p_for_sparsity[[s_lbl]]
      mask_seed <- run_id * 1000L + k
      masked <- apply_mask(cnts_base, p_sp, mask_seed)
      actual_sparsity <- round(sum(masked == 0L) / length(masked), 4)
      masked_cache[[s_lbl]] <- list(cnts = masked, p_sparsity_mask = p_sp,
                                     mask_seed = mask_seed, actual_sparsity = actual_sparsity)
    }

    # ----------------------------------------------------------------
    # Save one file per grid row sharing this fit
    # ----------------------------------------------------------------
    for (j in seq_len(nrow(rows))) {
      row_j    <- rows[j, ]
      out_path <- file.path(OUT_DIR,
        sprintf("symsim_run_%05d.rds", as.integer(row_j$run_id)))

      if (file.exists(out_path)) next

      mc   <- masked_cache[[as.character(row_j$sparsity)]]
      cnts <- mc$cnts

      cell_meta_out <- data.frame(
        run_id            = as.integer(row_j$run_id),
        cell_id           = paste0("cell_", seq_len(n_total)),
        true_group        = true_group,
        batch_id          = batch_id,
        achieved_sparsity = mc$actual_sparsity,
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
        sparsity_mask_p = mc$p_sparsity_mask,
        combined_mask_p = mc$p_sparsity_mask,
        mask_seed       = mc$mask_seed,
        actual_sparsity = mc$actual_sparsity,
        fit_key         = fk
      )

      saveRDS(
        object   = list(counts=cnts, cell_meta=cell_meta_out,
                        run_params=run_params),
        file     = out_path,
        compress = TRUE
      )
    }

    cat(sprintf("[fit %s] DONE  cells=%d  groups=%d  files=%d  combos=%d  %s\n",
      fk, n_total, length(unique(true_group)),
      nrow(rows), length(sparsity_labels_here), format(Sys.time(), "%H:%M:%S")))

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
