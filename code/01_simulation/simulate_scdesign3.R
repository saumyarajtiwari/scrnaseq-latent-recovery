# =============================================================================
# simulate_scdesign3.R  (REVISED — corrects a Step 1.7-discovered defect)
#
# CORRECTION NOTE:
#   The original version (preserved as simulate_scdesign3_v1_INVALID_sparsity_inert.R)
#   generated ONE count matrix per fit_key and saved it unchanged across all
#   5 sparsity_label values and all 3 dropout values sharing that fit_key.
#   Result: sparsity_label and dropout had ZERO effect on simulated data,
#   despite being logged in run_params as if they did. Discovered during
#   Step 1.7 output validation (byte-identical matrices across sparsity
#   labels within every one of 20 randomly sampled fit_key groups).
#
#   FIX: family_use is fixed at "nb" for all fits (faster and more reliable
#   than zinb, which was tested and found impractically slow — ~29 min/fit
#   vs ~5 min/fit for nb — and produced a broken all-NA zero-inflation
#   matrix on this scDesign3 version). dropout and sparsity_label are now
#   implemented as explicit, calibrated post-hoc stochastic zero-masking
#   applied to the generated count matrix:
#     - dropout_pi (0 / 0.10 / 0.40)  from scdesign3_calib_dropout.csv
#     - sparsity mask_p (0/.15/.35/.55/.75) from scdesign3_calib_sparsity.csv
#     - combined via: p_combined = 1 - (1 - dropout_pi) * (1 - mask_p)
#   Verified: strict monotonic sparsity ordering, zero degenerate all-zero
#   cells even at the harshest combined setting, tested at both extremes
#   of depth/separability (depth=500/separability=null and
#   depth=10000/separability=high).
#
#   Because dropout no longer requires a separate model fit, fit_key drops
#   it: unique fits fall from 244 to 82.
#
# Strategy: fit-once-mask-many
#   82 unique scDesign3 fits (depth x separability x n_cells x batch)
#   Each fit produces its full set of output files via masking + label
#   variation (dropout x sparsity x gene_strategy x clipping, where
#   gene_strategy and clipping remain preprocessing-stage labels only)
#
# Reference: PBMC 3k, top 2000 HVGs, annotated — data/simulated/pbmc3k_annotated.rds
# Output:    data/simulated/scdesign3/scdesign3_run_NNNNN.rds
# Log:       logs/scdesign3_run.log
#
# Estimated runtime : ~3-9 hours at 2 workers (real per-fit timing: 268-777 sec)
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

sparsity_calib <- read.csv("data/simulated/scdesign3_calib_sparsity.csv", stringsAsFactors=FALSE)
dropout_calib  <- read.csv("data/simulated/scdesign3_calib_dropout.csv", stringsAsFactors=FALSE)

mask_p_for_sparsity <- setNames(sparsity_calib$mask_p, as.character(sparsity_calib$sparsity_label))
pi_for_dropout       <- setNames(dropout_calib$dropout_pi, dropout_calib$dropout_label)

N_WORKERS <- 2L
OUT_DIR   <- "data/simulated/scdesign3"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

cat("=== scDesign3 Simulation (REVISED) ===\n")
cat(sprintf("Total runs       : %d\n",   nrow(param_grid)))
cat(sprintf("Parallel workers : %d\n",   N_WORKERS))
cat(sprintf("Output directory : %s\n",   OUT_DIR))
cat(sprintf("Started          : %s\n\n", format(Sys.time())))

# =============================================================================
# LOAD ANNOTATED REFERENCE
# =============================================================================

cat("Loading annotated reference...\n")
ref_full <- readRDS("data/simulated/pbmc3k_annotated.rds")
colnames(ref_full) <- paste0("cell_", seq_len(ncol(ref_full)))
cat(sprintf("Reference: %d genes x %d cells, %d cell types\n\n",
            nrow(ref_full), ncol(ref_full),
            length(unique(ref_full$cell_type))))

cat("Materializing counts into RAM (DelayedMatrix -> dgCMatrix)...\n")
counts(ref_full) <- as(counts(ref_full), "dgCMatrix")
gc()
cat("Materialization complete.\n\n")

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
# IDENTIFY UNIQUE FIT CONFIGURATIONS (dropout removed — now post-hoc)
# =============================================================================

fit_cols <- c("depth", "separability", "n_cells", "batch")
fit_key <- function(row) paste(row[fit_cols], collapse="_")

param_grid$fit_key <- apply(param_grid, 1, fit_key)
unique_fits        <- unique(param_grid$fit_key)
cat(sprintf("Unique scDesign3 fits: %d\n\n", length(unique_fits)))

# =============================================================================
# PER-FIT FUNCTION
# =============================================================================

run_one_fit <- function(fk) {

  rows <- param_grid[param_grid$fit_key == fk, ]
  row1 <- rows[1, ]

  out_paths <- file.path(OUT_DIR, sprintf("scdesign3_run_%05d.rds", rows$run_id))
  if (all(file.exists(out_paths))) {
    cat(sprintf("[fit %s] SKIP — all %d files exist\n", fk, nrow(rows)))
    return(invisible(NULL))
  }

  tryCatch({

    # 1. SUBSET REFERENCE BY SEPARABILITY
    sep_p      <- scdesign3_separability[[row1$separability]]
    cell_types <- sep_p$cell_types
    mu_formula <- sep_p$mu_formula

    keep_cells <- ref_full$cell_type %in% cell_types
    ref_sub    <- ref_full[, keep_cells]

    # 2. SCALE COUNTS BY DEPTH MULTIPLIER
    depth_p    <- scdesign3_depth[[as.character(row1$depth)]]
    mult       <- depth_p$lib_size_multiplier
    cnt        <- counts(ref_sub)
    cnt_scaled <- round(cnt * mult)
    counts(ref_sub) <- as(cnt_scaled, "dgCMatrix")

    # 3. SUBSAMPLE TO TARGET n_cells
    target_n <- as.integer(row1$n_cells)
    avail_n  <- ncol(ref_sub)
    actual_n <- min(target_n, avail_n)

    set.seed(as.integer(row1$run_id))
    keep_idx <- sample(avail_n, actual_n)
    ref_fit  <- ref_sub[, keep_idx]

    # 4. PREPARE FOR scDesign3
    colnames(ref_fit) <- paste0("cell_", seq_len(ncol(ref_fit)))
    ref_fit$cell_type <- factor(ref_fit$cell_type)
    counts(ref_fit)   <- as.matrix(counts(ref_fit))

    # 5. BATCH FORMULA
    if (row1$batch == "none") {
      final_formula <- mu_formula
    } else {
      if (row1$batch == "simple") {
        batch_vec <- rep(paste0("batch", 1:2),
                          c(floor(actual_n/2), actual_n - floor(actual_n/2)))
      } else {
        b1 <- round(actual_n * 0.50); b2 <- round(actual_n * 0.30); b3 <- actual_n - b1 - b2
        batch_vec <- rep(paste0("batch", 1:3), c(b1, b2, b3))
      }
      ref_fit$batch <- factor(batch_vec[seq_len(actual_n)])
      final_formula <- if (mu_formula == "1") "batch" else paste0(mu_formula, " + batch")
    }

    # 6. FIT scDesign3 (family_use ALWAYS "nb" — see correction note)
    other_cov <- if (row1$batch != "none") "batch" else NULL

    fit <- scdesign3(
      sce = ref_fit, assay_use = "counts", celltype = "cell_type",
      pseudotime = NULL, spatial = NULL, other_covariates = other_cov,
      mu_formula = final_formula, sigma_formula = "1", family_use = "nb",
      n_cores = 1L, usebam = FALSE, corr_formula = "1", copula = "gaussian",
      DT = TRUE, pseudo_obs = FALSE, return_model = FALSE
    )

    if (is.null(fit) || is.null(fit$new_count)) stop("scdesign3() returned NULL or empty new_count")

    cnts_base <- as(fit$new_count, "dgCMatrix")

    new_cov_df <- fit$new_covariate
    true_group <- if ("cell_type" %in% names(new_cov_df)) as.character(new_cov_df$cell_type) else rep(cell_types[1], ncol(cnts_base))
    batch_id   <- if ("batch" %in% names(new_cov_df)) as.integer(sub("batch", "", as.character(new_cov_df$batch))) else rep(1L, ncol(cnts_base))

    # 7. FOR EACH (dropout, sparsity_label) PAIR: APPLY COMBINED MASKING ONCE
    combos <- unique(rows[, c("dropout", "sparsity")])
    masked_cache <- list()

    for (k in seq_len(nrow(combos))) {
      d_lvl  <- combos$dropout[k]
      s_lbl  <- as.character(combos$sparsity[k])
      p_drop <- pi_for_dropout[[d_lvl]]
      p_sp   <- mask_p_for_sparsity[[s_lbl]]
      p_comb <- 1 - (1 - p_drop) * (1 - p_sp)

      mask_seed <- as.integer(row1$run_id) * 1000L + k
      masked <- apply_mask(cnts_base, p_comb, mask_seed)
      actual_sparsity <- round(1 - Matrix::nnzero(masked) / length(masked), 4)

      masked_cache[[paste(d_lvl, s_lbl, sep="__")]] <- list(
        cnts = masked, p_combined = p_comb, p_dropout = p_drop,
        p_sparsity_mask = p_sp, actual_sparsity = actual_sparsity,
        mask_seed = mask_seed
      )
    }

    # 8. SAVE ONE FILE PER GRID ROW
    for (j in seq_len(nrow(rows))) {
      row_j    <- rows[j, ]
      out_path <- file.path(OUT_DIR, sprintf("scdesign3_run_%05d.rds", row_j$run_id))
      if (file.exists(out_path)) next

      key  <- paste(row_j$dropout, as.character(row_j$sparsity), sep="__")
      mc   <- masked_cache[[key]]
      cnts <- mc$cnts

      cell_meta <- data.frame(
        run_id            = row_j$run_id,
        cell_id           = colnames(cnts),
        true_group        = true_group,
        batch_id          = batch_id,
        achieved_sparsity = mc$actual_sparsity,
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
        family_use          = "nb",
        depth_multiplier    = mult,
        dropout_pi          = mc$p_dropout,
        sparsity_mask_p     = mc$p_sparsity_mask,
        combined_mask_p     = mc$p_combined,
        mask_seed           = mc$mask_seed,
        actual_sparsity     = mc$actual_sparsity,
        fit_key             = fk
      )

      saveRDS(
        object   = list(counts=cnts, cell_meta=cell_meta, run_params=run_params),
        file     = out_path,
        compress = TRUE
      )
    }

    cat(sprintf("[fit %s] DONE  cells=%d  files=%d  combos=%d  %s\n",
                fk, ncol(cnts_base), nrow(rows), nrow(combos),
                format(Sys.time(), "%H:%M:%S")))

  }, error = function(e) {
    cat(sprintf("[fit %s] ERROR  %s\n", fk, conditionMessage(e)))
  })
}

# =============================================================================
# RUN ALL FITS (parallel, 2 workers — nb confirmed fork-safe in original run)
# =============================================================================

cat("Starting generation...\n\n")
invisible(mclapply(unique_fits, run_one_fit, mc.cores = N_WORKERS))

cat(sprintf("\n=== Complete: %s ===\n", format(Sys.time())))
