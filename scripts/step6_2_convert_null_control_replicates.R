# Convert the 30 null-control replicate files (rep2, rep3; Step 1.6) from their
# raw list format (counts, cell_meta, run_params) into standardized
# SingleCellExperiment objects matching Step 1.8's schema, so the existing,
# unmodified Step 3 preprocessing scripts can discover and process them via
# their normal data/simulated/sce/{simulator}/ scan.
#
# Scope note: true_group_means / ground_truth_source metadata (populated for
# the main grid in Steps 1.9-1.12) is deliberately NOT reconstructed here.
# Ground-truth extraction was never run against these 30 replicates, and
# nothing currently planned (Step 6.3, or any completed step) requires it —
# Step 4's subspace-recovery metrics already exclude all null-controls
# entirely (rank-0 true subspace is undefined). Documented here rather than
# silently fabricated.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

src_dir  <- "data/simulated/null_control"
out_base <- "data/simulated/sce"

simulators <- c("splatter", "scdesign3", "symsim")
converted <- 0
skipped <- 0
errors <- list()

for (sim in simulators) {
  sim_dir <- file.path(src_dir, sim)
  files <- list.files(sim_dir, pattern = "rep[23]\\.rds$", full.names = TRUE)
  cat(sprintf("\n=== %s: %d replicate files found ===\n", sim, length(files)))

  for (fp in files) {
    res <- tryCatch({
      raw <- readRDS(fp)
      cm  <- raw$cell_meta
      rp  <- raw$run_params

      run_id <- unique(cm$run_id)
      rep_n  <- unique(cm$replicate)
      if (length(run_id) != 1 || length(rep_n) != 1) {
        stop("expected single run_id/replicate value per file, found multiple")
      }

      out_file <- file.path(out_base, sim,
                             sprintf("%s_sce_run_%05d_rep%d.rds", sim, run_id, rep_n))

      if (file.exists(out_file)) {
        return(list(status = "skipped_existing", out_file = out_file))
      }

      col_data <- DataFrame(
        run_id = cm$run_id,
        cell_id = cm$cell_id,
        true_group = cm$true_group,
        batch_id = cm$batch_id,
        achieved_sparsity = cm$achieved_sparsity
      )

      sce <- SingleCellExperiment(
        assays = list(counts = raw$counts),
        colData = col_data
      )
      metadata(sce)$run_params <- rp
      metadata(sce)$simulator <- sim
      metadata(sce)$conversion_note <- paste0(
        "Converted from raw null-control replicate format (",
        basename(fp), ") on ", Sys.Date(),
        " for Step 6.3. true_group_means/ground_truth_source deliberately ",
        "omitted (never extracted for these replicates; not required by ",
        "any currently planned step)."
      )

      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
      saveRDS(sce, tmp_file)
      file.rename(tmp_file, out_file)

      list(status = "converted", out_file = out_file, n_cells = ncol(sce), n_genes = nrow(sce))
    }, error = function(e) {
      list(status = "error", msg = conditionMessage(e), file = fp)
    })

    if (res$status == "converted") {
      converted <- converted + 1
      cat(sprintf("  OK: %s -> %s (%d genes x %d cells)\n",
                   basename(fp), basename(res$out_file), res$n_genes, res$n_cells))
    } else if (res$status == "skipped_existing") {
      skipped <- skipped + 1
      cat(sprintf("  SKIP (exists): %s\n", basename(fp)))
    } else {
      errors[[length(errors) + 1]] <- res
      cat(sprintf("  ERROR: %s -> %s\n", basename(fp), res$msg))
    }
  }
}

cat(sprintf("\n=== Summary: %d converted, %d skipped, %d errors ===\n",
            converted, skipped, length(errors)))
if (length(errors) > 0) {
  for (e in errors) cat("  ", e$file, ":", e$msg, "\n")
}
