# Step 1.6 Phase B (scDesign3 only) — isolated per original architecture pattern
# REVISED — see simulate_scdesign3.R correction note for full context.
# Fix: sparsity_label now applies calibrated post-hoc masking (dropout is
# fixed at "none"/pi=0 for null-control by design, so only the sparsity
# ladder applies here).
suppressPackageStartupMessages({ library(scDesign3); library(SingleCellExperiment); library(Matrix) })
source("R/01_simulation/param_dict.R")
manifest <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
rows <- manifest[manifest$is_new == TRUE & manifest$simulator == "scdesign3", ]
cat("=== Null-control Phase B: scDesign3 (REVISED) ===\n")

sparsity_calib <- read.csv("data/simulated/scdesign3_calib_sparsity.csv", stringsAsFactors=FALSE)
mask_p_for_sparsity <- setNames(sparsity_calib$mask_p, as.character(sparsity_calib$sparsity_label))

apply_mask <- function(mat, p, seed) {
  if (p <= 0) return(mat)
  set.seed(seed)
  mat_t <- as(mat, "TsparseMatrix")
  n_nz <- length(mat_t@x)
  keep <- rbinom(n_nz, 1, 1 - p) == 1
  mat_t@x[!keep] <- 0
  drop0(as(mat_t, "CsparseMatrix"))
}

ref_full <- readRDS("data/simulated/pbmc3k_annotated.rds")
counts(ref_full) <- as(counts(ref_full), "dgCMatrix")
for (rep in unique(rows$replicate)) {
  rep_rows <- rows[rows$replicate == rep, ]
  if (all(file.exists(rep_rows$file_path))) { cat("[scdesign3] rep", rep, "SKIP\n"); next }
  seed <- rep_rows$seed[1]
  t0 <- Sys.time()
  sep_p      <- scdesign3_separability[["null"]]
  keep_cells <- ref_full$cell_type %in% sep_p$cell_types
  ref_sub    <- ref_full[, keep_cells]
  depth_p    <- scdesign3_depth[["2000"]]
  cnt_scaled <- round(counts(ref_sub) * depth_p$lib_size_multiplier)
  counts(ref_sub) <- as(cnt_scaled, "dgCMatrix")
  set.seed(seed)
  keep_idx <- sample(ncol(ref_sub), min(1000L, ncol(ref_sub)))
  ref_fit  <- ref_sub[, keep_idx]
  colnames(ref_fit) <- paste0("cell_", seq_len(ncol(ref_fit)))
  ref_fit$cell_type <- factor(ref_fit$cell_type)
  counts(ref_fit)   <- as.matrix(counts(ref_fit))
  fit <- scdesign3(
    sce = ref_fit, assay_use = "counts", celltype = "cell_type",
    pseudotime = NULL, spatial = NULL, other_covariates = NULL,
    mu_formula = sep_p$mu_formula, sigma_formula = "1",
    family_use = "nb", n_cores = 1L, usebam = FALSE,
    corr_formula = "1", copula = "gaussian", DT = TRUE,
    pseudo_obs = FALSE, return_model = FALSE
  )
  if (is.null(fit) || is.null(fit$new_count)) stop("scdesign3() returned NULL/empty new_count")
  cnts_base <- as(fit$new_count, "dgCMatrix")
  new_cov_df <- fit$new_covariate
  true_group <- if ("cell_type" %in% names(new_cov_df)) as.character(new_cov_df$cell_type) else rep("CD4_T", ncol(cnts_base))

  for (j in seq_len(nrow(rep_rows))) {
    rj <- rep_rows[j, ]
    if (file.exists(rj$file_path)) next

    s_lbl <- as.character(rj$sparsity)
    p_sp  <- mask_p_for_sparsity[[s_lbl]]
    mask_seed <- as.integer(rj$base_run_id) * 1000L + as.integer(rj$replicate) * 10L + j
    cnts <- apply_mask(cnts_base, p_sp, mask_seed)
    actual_sparsity <- round(1 - Matrix::nnzero(cnts) / length(cnts), 4)

    cell_meta <- data.frame(
      run_id = rj$base_run_id, replicate = rj$replicate,
      cell_id = colnames(cnts), true_group = true_group, batch_id = 1L,
      achieved_sparsity = actual_sparsity, stringsAsFactors = FALSE
    )
    run_params <- list(
      run_id = rj$base_run_id, replicate = rj$replicate, seed = seed,
      sparsity_label = s_lbl, depth_label = "2000",
      dropout = "none", separability = "null", n_cells_target = 1000L,
      n_cells_actual = ncol(cnts), n_groups = length(unique(true_group)),
      batch = "none", gene_strategy = "all", clipping = "none",
      is_null_control = TRUE,
      family_use = "nb", depth_multiplier = depth_p$lib_size_multiplier,
      dropout_pi = 0, sparsity_mask_p = p_sp, combined_mask_p = p_sp,
      mask_seed = mask_seed, actual_sparsity = actual_sparsity
    )
    dir.create(dirname(rj$file_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
            rj$file_path, compress = TRUE)
  }
  cat(sprintf("[scdesign3] rep=%d DONE cells=%d (%.1fs)\n",
              rep, ncol(cnts_base), as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("=== scDesign3 Phase B complete: %s ===\n", format(Sys.time())))
