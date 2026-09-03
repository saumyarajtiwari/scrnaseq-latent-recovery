# Step 1.6 Phase B (SymSim only) — isolated per original architecture pattern
# REVISED — see simulate_symsim.R correction note for full context.
# Fix: sparsity_label now applies calibrated post-hoc masking (dropout is
# fixed at "none" for null-control by design, unaffected by this change).
suppressPackageStartupMessages({ library(SymSim); library(Matrix) })
source("R/01_simulation/param_dict.R")
manifest <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
sym_rows <- manifest[manifest$is_new == TRUE & manifest$simulator == "symsim", ]
cat("=== Null-control Phase B: SymSim (REVISED) ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

sparsity_calib <- read.csv("data/simulated/symsim_calib_sparsity.csv", stringsAsFactors=FALSE)
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

for (rep in unique(sym_rows$replicate)) {
  rep_rows <- sym_rows[sym_rows$replicate == rep, ]
  if (all(file.exists(rep_rows$file_path))) { cat("[symsim] rep", rep, "SKIP\n"); next }
  seed <- as.integer(rep_rows$seed[1])
  t0 <- Sys.time()
  cat(sprintf("[symsim] rep=%d starting SimulateTrueCounts at %s...\n", rep, format(t0)))
  depth_p   <- symsim_depth[["none"]][["2000"]]
  dropout_p <- symsim_dropout[["none"]]
  tr <- SimulateTrueCounts(
    ncells_total = 1000L, min_popsize = 1000L, i_minpop = 1L,
    ngenes = 2000L, evf_type = "one.population", phyla = NULL,
    randseed = seed, Sigma = SYMSIM_SIGMA, n_de_evf = 8L, bimod = 0, vary = "s"
  )
  cat(sprintf("  SimulateTrueCounts done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
  true_counts  <- tr[["counts"]]
  cell_meta_tr <- tr[["cell_meta"]]
  set.seed(seed)
  obs <- True2ObservedCounts(
    true_counts = true_counts, meta_cell = cell_meta_tr, protocol = "UMI",
    alpha_mean = dropout_p$alpha_mean, alpha_sd = dropout_p$alpha_sd,
    gene_len = rep(1000L, 2000L),
    depth_mean = depth_p$depth_mean, depth_sd = depth_p$depth_sd
  )
  cnts_base <- as(obs[[1]], "dgCMatrix")
  true_group <- as.character(cell_meta_tr$pop)

  for (j in seq_len(nrow(rep_rows))) {
    rj <- rep_rows[j, ]
    if (file.exists(rj$file_path)) next

    s_lbl <- as.character(rj$sparsity)
    p_sp  <- mask_p_for_sparsity[[s_lbl]]
    mask_seed <- as.integer(rj$base_run_id) * 1000L + as.integer(rj$replicate) * 10L + j
    cnts <- apply_mask(cnts_base, p_sp, mask_seed)
    actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)

    cell_meta <- data.frame(
      run_id = rj$base_run_id, replicate = rj$replicate,
      cell_id = paste0("cell_", seq_len(ncol(cnts))),
      true_group = true_group, batch_id = 1L,
      achieved_sparsity = actual_sparsity, stringsAsFactors = FALSE
    )
    run_params <- list(
      run_id = rj$base_run_id, replicate = rj$replicate, seed = seed,
      sparsity_label = s_lbl, depth_label = "2000",
      dropout = "none", separability = "null", n_cells = 1000L,
      n_cells_actual = ncol(cnts), n_groups = length(unique(true_group)),
      batch = "none", gene_strategy = "all", clipping = "none", is_null_control = TRUE,
      sigma = SYMSIM_SIGMA, n_de_evf = 8L, alpha_mean = dropout_p$alpha_mean,
      depth_mean = depth_p$depth_mean,
      sparsity_mask_p = p_sp, combined_mask_p = p_sp, mask_seed = mask_seed,
      actual_sparsity = actual_sparsity
    )
    dir.create(dirname(rj$file_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
            rj$file_path, compress = TRUE)
  }
  cat(sprintf("[symsim] rep=%d DONE cells=%d groups=%d (%.1fs total)\n",
              rep, ncol(cnts_base), length(unique(true_group)),
              as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("\n=== SymSim Phase B complete: %s ===\n", format(Sys.time())))
