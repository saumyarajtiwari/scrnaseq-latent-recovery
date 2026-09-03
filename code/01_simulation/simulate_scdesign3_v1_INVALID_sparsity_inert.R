# =============================================================================
# simulate_scdesign3.R
#
# Runs all 10,940 scDesign3 simulations from param_grid.csv.
# One .rds output per run: sparse count matrix + cell metadata + run parameters.
#
# Strategy: fit-once-simulate-many
#   244 unique scDesign3 fits (depth x dropout x separability x n_cells x batch)
#   Each fit produces 45 output files (main) or 5 (null) via metadata variation
#   sparsity labels are ordinal identifiers — do not change the count matrix
#   gene_strategy and clipping are preprocessing-stage labels only
#
# Reference: PBMC 3k, top 2000 HVGs, annotated — data/simulated/pbmc3k_annotated.rds
# Output:    data/simulated/scdesign3/scdesign3_run_NNNNN.rds
# Log:       logs/scdesign3_run.log
#
# Estimated runtime : 15-25 hours at 2 workers
# =============================================================================

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

N_WORKERS <- 2L
OUT_DIR   <- "data/simulated/scdesign3"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== scDesign3 Simulation ===\n")
cat(sprintf("Total runs       : %d\n",   nrow(param_grid)))
cat(sprintf("Parallel workers : %d\n",   N_WORKERS))
cat(sprintf("Output directory : %s\n",   OUT_DIR))
cat(sprintf("Started          : %s\n\n", format(Sys.time())))

# =============================================================================
# LOAD ANNOTATED REFERENCE
# =============================================================================

cat("Loading annotated reference...\n")
ref_full <- readRDS("data/simulated/pbmc3k_annotated.rds")
cat(sprintf("Reference: %d genes x %d cells, %d cell types\n\n",
            nrow(ref_full), ncol(ref_full),
            length(unique(ref_full$cell_type))))

cat("Materializing counts into RAM (DelayedMatrix -> dgCMatrix)...\n")
suppressPackageStartupMessages(library(SingleCellExperiment))
counts(ref_full) <- as(counts(ref_full), "dgCMatrix")
gc()
cat("Materialization complete.\n\n")

# =============================================================================
# IDENTIFY UNIQUE FIT CONFIGURATIONS
# =============================================================================

# Columns that determine a unique scDesign3 fit
fit_cols <- c("depth", "dropout", "separability", "n_cells", "batch")

fit_key <- function(row) {
  paste(row[fit_cols], collapse="_")
}

param_grid$fit_key <- apply(param_grid, 1, fit_key)
unique_fits        <- unique(param_grid$fit_key)
cat(sprintf("Unique scDesign3 fits: %d\n\n", length(unique_fits)))

# =============================================================================
# PER-FIT FUNCTION
# =============================================================================

run_one_fit <- function(fk) {

  rows   <- param_grid[param_grid$fit_key == fk, ]
  row1   <- rows[1, ]

  # Check if all output files for this fit already exist (checkpointing)
  out_paths <- file.path(OUT_DIR,
    sprintf("scdesign3_run_%05d.rds", rows$run_id))
  if (all(file.exists(out_paths))) {
    cat(sprintf("[fit %s] SKIP — all %d files exist\n", fk, nrow(rows)))
    return(invisible(NULL))
  }

  tryCatch({

    # ------------------------------------------------------------------
    # 1. SUBSET REFERENCE BY SEPARABILITY
    # ------------------------------------------------------------------
    sep_p      <- scdesign3_separability[[row1$separability]]
    cell_types <- sep_p$cell_types
    mu_formula <- sep_p$mu_formula

    keep_cells <- ref_full$cell_type %in% cell_types
    ref_sub    <- ref_full[, keep_cells]

    # ------------------------------------------------------------------
    # 2. SCALE COUNTS BY DEPTH MULTIPLIER
    # ------------------------------------------------------------------
    depth_p    <- scdesign3_depth[[as.character(row1$depth)]]
    mult       <- depth_p$lib_size_multiplier
    cnt        <- counts(ref_sub)
    cnt_scaled <- round(cnt * mult)
    cnt_scaled <- as(cnt_scaled, "dgCMatrix")
    counts(ref_sub) <- cnt_scaled

    # ------------------------------------------------------------------
    # 3. SUBSAMPLE TO TARGET n_cells
    # ------------------------------------------------------------------
    target_n   <- as.integer(row1$n_cells)
    avail_n    <- ncol(ref_sub)
    actual_n   <- min(target_n, avail_n)

    set.seed(as.integer(row1$run_id))
    keep_idx   <- sample(avail_n, actual_n)
    ref_fit    <- ref_sub[, keep_idx]

    # ------------------------------------------------------------------
    # 4. PREPARE FOR scDesign3
    # ------------------------------------------------------------------
    colnames(ref_fit) <- paste0("cell_", seq_len(ncol(ref_fit)))
    ref_fit$cell_type <- factor(ref_fit$cell_type)
    counts(ref_fit)   <- as.matrix(counts(ref_fit))

    # ------------------------------------------------------------------
    # 5. DROPOUT / FAMILY
    # ------------------------------------------------------------------
    dropout_p  <- scdesign3_dropout[[row1$dropout]]
    family_use <- dropout_p$family_use

    # ------------------------------------------------------------------
    # 6. BATCH FORMULA
    # ------------------------------------------------------------------
    if (row1$batch == "none") {
      final_formula <- mu_formula
      n_batches     <- 1L
    } else {
      # Add batch covariate to reference
      if (row1$batch == "simple") {
        n_batches  <- 2L
        batch_vec  <- rep(paste0("batch", 1:2),
                          c(floor(actual_n/2), actual_n - floor(actual_n/2)))
      } else {
        n_batches  <- 3L
        b1 <- round(actual_n * 0.50)
        b2 <- round(actual_n * 0.30)
        b3 <- actual_n - b1 - b2
        batch_vec  <- rep(paste0("batch", 1:3), c(b1, b2, b3))
      }
      ref_fit$batch <- factor(batch_vec[seq_len(actual_n)])
      if (mu_formula == "1") {
        final_formula <- "batch"
      } else {
        final_formula <- paste0(mu_formula, " + batch")
      }
    }

    # ------------------------------------------------------------------
    # 7. FIT scDesign3
    # ------------------------------------------------------------------
    other_cov <- if (row1$batch != "none") "batch" else NULL

    fit <- scdesign3(
      sce              = ref_fit,
      assay_use        = "counts",
      celltype         = "cell_type",
      pseudotime       = NULL,
      spatial          = NULL,
      other_covariates = other_cov,
      mu_formula       = final_formula,
      sigma_formula    = "1",
      family_use       = family_use,
      n_cores          = 1L,
      usebam           = FALSE,
      corr_formula     = "1",
      copula           = "gaussian",
      DT               = TRUE,
      pseudo_obs       = FALSE,
      return_model     = FALSE
    )

    if (is.null(fit) || is.null(fit$new_count)) {
      stop("scdesign3() returned NULL or empty new_count")
    }

    # Generated count matrix
    cnts           <- as(fit$new_count, "dgCMatrix")
    actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)

    # Cell metadata from new_covariate output
    new_cov_df     <- fit$new_covariate
    true_group     <- if ("cell_type" %in% names(new_cov_df)) {
      as.character(new_cov_df$cell_type)
    } else {
      rep(cell_types[1], ncol(cnts))
    }
    batch_id       <- if ("batch" %in% names(new_cov_df)) {
      as.integer(sub("batch", "", as.character(new_cov_df$batch)))
    } else {
      rep(1L, ncol(cnts))
    }

    # ------------------------------------------------------------------
    # 8. SAVE ONE FILE PER GRID ROW SHARING THIS FIT
    # ------------------------------------------------------------------
    for (j in seq_len(nrow(rows))) {
      row_j    <- rows[j, ]
      out_path <- file.path(OUT_DIR,
        sprintf("scdesign3_run_%05d.rds", row_j$run_id))

      if (file.exists(out_path)) next

      cell_meta <- data.frame(
        run_id            = row_j$run_id,
        cell_id           = colnames(cnts),
        true_group        = true_group,
        batch_id          = batch_id,
        achieved_sparsity = actual_sparsity,
        stringsAsFactors  = FALSE
      )

      run_params <- list(
        run_id              = row_j$run_id,
        sparsity_label      = as.character(row_j$sparsity),
        depth_label         = as.character(row_j$depth),
        dropout             = row_j$dropout,
        separability        = row_j$separability,
        n_cells_target      = as.integer(row_j$n_cells),
        n_cells_actual      = ncol(cnts),
        n_groups            = length(unique(true_group)),
        batch               = row_j$batch,
        gene_strategy       = row_j$gene_strategy,
        clipping            = row_j$clipping,
        is_null_control     = as.logical(row_j$is_null_control),
        family_use          = family_use,
        depth_multiplier    = mult,
        actual_sparsity     = actual_sparsity,
        fit_key             = fk
      )

      saveRDS(
        object   = list(counts=cnts, cell_meta=cell_meta, run_params=run_params),
        file     = out_path,
        compress = TRUE
      )
    }

    cat(sprintf("[fit %s] DONE  sparsity=%.4f  cells=%d  files=%d  %s\n",
                fk, actual_sparsity, ncol(cnts), nrow(rows),
                format(Sys.time(), "%H:%M:%S")))

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

n_done   <- length(list.files(OUT_DIR, pattern="^scdesign3_run_.*\\.rds$"))
n_target <- nrow(param_grid)

cat(sprintf("\n=== Run complete: %s ===\n", format(Sys.time())))
cat(sprintf("Completed : %d / %d\n", n_done, n_target))
if (n_done < n_target) {
  cat(sprintf("Failed    : %d  — re-run to retry\n", n_target - n_done))
}
