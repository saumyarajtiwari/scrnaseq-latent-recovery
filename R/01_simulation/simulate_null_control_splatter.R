# Step 1.6 Phase B (Splatter only) — isolated per original architecture pattern
suppressPackageStartupMessages({ library(splatter); library(Matrix) })
source("R/01_simulation/param_dict.R")

manifest <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
rows <- manifest[manifest$is_new == TRUE & manifest$simulator == "splatter", ]

cat("=== Null-control Phase B: Splatter ===\n")
for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  if (file.exists(r$file_path)) { cat("[splatter]", r$file_path, "SKIP\n"); next }

  sparsity_p <- splatter_sparsity[[sprintf("%.2f", as.numeric(r$sparsity))]]
  depth_p    <- splatter_depth[["2000"]]
  t0 <- Sys.time()

  sim <- splatSimulate(
    nGenes = 10000L, batchCells = 1000L,
    bcv.common = sparsity_p$bcv.common,
    lib.loc = depth_p$lib.loc, lib.scale = depth_p$lib.scale,
    dropout.type = "none", method = "single",
    batch.facLoc = 0.1, batch.facScale = 0.1,
    verbose = FALSE, seed = r$seed
  )
  cnts <- as(counts(sim), "dgCMatrix")
  actual_sparsity <- round(sum(cnts == 0L) / length(cnts), 4)

  cell_meta <- data.frame(
    run_id = r$base_run_id, replicate = r$replicate,
    cell_id = colnames(sim), true_group = "Group1", batch_id = 1L,
    achieved_sparsity = actual_sparsity, stringsAsFactors = FALSE
  )
  run_params <- list(
    run_id = r$base_run_id, replicate = r$replicate, seed = r$seed,
    sparsity_label = as.character(r$sparsity), depth_label = "2000",
    dropout = "none", separability = "null", n_cells = 1000L, n_groups = 1L,
    batch = "none", gene_strategy = "all", clipping = "none", is_null_control = TRUE,
    bcv_common = sparsity_p$bcv.common, lib_loc = depth_p$lib.loc,
    actual_sparsity = actual_sparsity
  )
  dir.create(dirname(r$file_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(counts = cnts, cell_meta = cell_meta, run_params = run_params),
          r$file_path, compress = TRUE)
  cat(sprintf("[splatter] sparsity=%.2f rep=%d DONE actual=%.4f (%.1fs)\n",
              r$sparsity, r$replicate, actual_sparsity,
              as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("=== Splatter Phase B complete: %s ===\n", format(Sys.time())))
