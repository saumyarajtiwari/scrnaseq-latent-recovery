# Step 1.6 Phase B (SymSim only) — isolated per original architecture pattern
suppressPackageStartupMessages({ library(SymSim); library(Matrix) })
source("R/01_simulation/param_dict.R")

manifest <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
sym_rows <- manifest[manifest$is_new == TRUE & manifest$simulator == "symsim", ]

cat("=== Null-control Phase B: SymSim ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

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

  cnts <- as(obs[[1]], "dgCMatrix")
  actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)
  true_group <- as.character(cell_meta_tr$pop)

  for (j in seq_len(nrow(rep_rows))) {
    rj <- rep_rows[j, ]
    if (file.exists(rj$file_path)) next
    cell_meta <- data.frame(
      run_id = rj$base_run_id, replicate = rj$replicate,
      cell_id = paste0("cell_", seq_len(ncol(cnts))),
      true_group = true_group, batch_id = 1L,
      achieved_sparsity = actual_sparsity, stringsAsFactors = FALSE
    )
    run_params <- list(
      run_id = rj$base_run_id, replicate = rj$replicate, seed = seed,
      sparsity_label = as.character(rj$sparsity), depth_label = "2000",
      dropout = "none", separability = "null", n_cells = 1000L,
      n_cells_actual = ncol(cnts), n_groups = length(unique(true_group)),
      batch = "none", is_null_control = TRUE,
      sigma = SYMSIM_SIGMA, n_de_evf = 8L, alpha_mean = dropout_p$alpha_mean,
      depth_mean = depth_p$depth_mean, actual_sparsity = actual_sparsity
    )
    dir.create(dirname(rj$file_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
            rj$file_path, compress = TRUE)
  }
  cat(sprintf("[symsim] rep=%d DONE actual=%.4f cells=%d groups=%d (%.1fs total)\n",
              rep, actual_sparsity, ncol(cnts), length(unique(true_group)),
              as.numeric(Sys.time() - t0, units = "secs")))
}

cat(sprintf("\n=== SymSim Phase B complete: %s ===\n", format(Sys.time())))
