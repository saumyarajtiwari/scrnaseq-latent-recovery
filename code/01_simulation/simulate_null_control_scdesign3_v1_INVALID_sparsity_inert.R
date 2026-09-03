# Step 1.6 Phase B (scDesign3 only) — isolated per original architecture pattern
suppressPackageStartupMessages({ library(scDesign3); library(SingleCellExperiment); library(Matrix) })
source("R/01_simulation/param_dict.R")

manifest <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
rows <- manifest[manifest$is_new == TRUE & manifest$simulator == "scdesign3", ]

cat("=== Null-control Phase B: scDesign3 ===\n")
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

  dropout_p <- scdesign3_dropout[["none"]]

  fit <- scdesign3(
    sce = ref_fit, assay_use = "counts", celltype = "cell_type",
    pseudotime = NULL, spatial = NULL, other_covariates = NULL,
    mu_formula = sep_p$mu_formula, sigma_formula = "1",
    family_use = dropout_p$family_use, n_cores = 1L, usebam = FALSE,
    corr_formula = "1", copula = "gaussian", DT = TRUE,
    pseudo_obs = FALSE, return_model = FALSE
  )
  if (is.null(fit) || is.null(fit$new_count)) stop("scdesign3() returned NULL/empty new_count")

  cnts <- as(fit$new_count, "dgCMatrix")
  actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)
  new_cov_df <- fit$new_covariate
  true_group <- if ("cell_type" %in% names(new_cov_df)) as.character(new_cov_df$cell_type) else rep("CD4_T", ncol(cnts))

  for (j in seq_len(nrow(rep_rows))) {
    rj <- rep_rows[j, ]
    if (file.exists(rj$file_path)) next
    cell_meta <- data.frame(
      run_id = rj$base_run_id, replicate = rj$replicate,
      cell_id = colnames(cnts), true_group = true_group, batch_id = 1L,
      achieved_sparsity = actual_sparsity, stringsAsFactors = FALSE
    )
    run_params <- list(
      run_id = rj$base_run_id, replicate = rj$replicate, seed = seed,
      sparsity_label = as.character(rj$sparsity), depth_label = "2000",
      dropout = "none", separability = "null", n_cells_target = 1000L,
      n_cells_actual = ncol(cnts), n_groups = length(unique(true_group)),
      batch = "none", is_null_control = TRUE,
      family_use = dropout_p$family_use, depth_multiplier = depth_p$lib_size_multiplier,
      actual_sparsity = actual_sparsity
    )
    dir.create(dirname(rj$file_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
            rj$file_path, compress = TRUE)
  }
  cat(sprintf("[scdesign3] rep=%d DONE actual=%.4f cells=%d (%.1fs)\n",
              rep, actual_sparsity, ncol(cnts), as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("=== scDesign3 Phase B complete: %s ===\n", format(Sys.time())))
