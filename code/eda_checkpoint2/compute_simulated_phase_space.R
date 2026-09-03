# EDA Checkpoint 2.2 (part 1) — Compute achieved sparsity/depth for every
# simulated run, to build the phase-space background against which real
# datasets will be overlaid. Checkpoints every 2000 files processed.

suppressPackageStartupMessages(library(SingleCellExperiment))

files <- list.files("data/simulated/sce", pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
cat("Total simulated files found:", length(files), "\n")

checkpoint_path <- "results/eda_checkpoint2/simulated_phase_space_summary.csv"
dir.create(dirname(checkpoint_path), recursive = TRUE, showWarnings = FALSE)

results <- vector("list", length(files))
start_idx <- 1

if (file.exists(checkpoint_path)) {
  existing <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
  cat("Found existing checkpoint with", nrow(existing), "rows. Resuming after that point.\n")
  start_idx <- nrow(existing) + 1
  for (i in seq_len(nrow(existing))) results[[i]] <- existing[i, ]
}

for (i in start_idx:length(files)) {
  f <- files[i]
  simulator <- basename(dirname(f))

  sce <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(sce)) {
    cat("SKIPPED (unreadable):", f, "\n")
    next
  }

  lib_size <- Matrix::colSums(assay(sce, "counts"))
  rp <- metadata(sce)$run_params

  results[[i]] <- data.frame(
    run_id = rp$run_id,
    simulator = simulator,
    sparsity_label = rp$sparsity_label,
    depth_label = rp$depth_label,
    dropout = rp$dropout,
    separability = rp$separability,
    is_null_control = rp$is_null_control,
    n_cells = ncol(sce),
    mean_lib_size = mean(lib_size),
    achieved_sparsity_pct = 100 * unique(sce$achieved_sparsity)[1],
    stringsAsFactors = FALSE
  )

  if (i %% 2000 == 0 || i == length(files)) {
    cat("Processed", i, "of", length(files), "-- checkpointing...\n")
    combined <- do.call(rbind, results[1:i])
    write.csv(combined, checkpoint_path, row.names = FALSE)
  }
}

cat("Done. Final file:", checkpoint_path, "\n")
