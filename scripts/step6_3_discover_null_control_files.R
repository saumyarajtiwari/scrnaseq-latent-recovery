# Purpose-built file discovery for Step 6.3 (Phantom Clustering), scoped only
# to null-control files. Does NOT modify or re-run step3_10_manifest.R, whose
# filename regex is anchored to the main-grid pattern and cannot parse the
# new _rep2/_rep3 suffix introduced when closing the Step 1.6/Step 3 gap.

METHOD_MAP <- list(
  pca_raw = "rawpca", pca_libnorm = "libnormpca", pca_log = "logpca",
  pca_shiftedlog = "shiftedlogpca", pca_sctransform_v2 = "sctpca", pca_glmpca = "glmpca"
)

nc_grid <- read.csv("data/simulated/null_control_grid.csv", stringsAsFactors = FALSE)
cat("null_control_grid.csv rows (should be 45):", nrow(nc_grid), "\n")
cat("Columns:", paste(colnames(nc_grid), collapse=", "), "\n")

rows <- list()
for (method_dir in names(METHOD_MAP)) {
  methodtag <- METHOD_MAP[[method_dir]]
  for (sim in c("splatter", "scdesign3", "symsim")) {
    d <- file.path("data/processed", method_dir, "simulated", sim)
    files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)

    pat_rep1 <- paste0("^", sim, "_sce_run_([0-9]{5})_", methodtag, "\\.rds$")
    pat_repN <- paste0("^", sim, "_sce_run_([0-9]{5})_rep([23])_", methodtag, "\\.rds$")

    for (f in files) {
      bn <- basename(f)
      m1 <- regmatches(bn, regexec(pat_rep1, bn))[[1]]
      mN <- regmatches(bn, regexec(pat_repN, bn))[[1]]
      if (length(m1) == 2) {
        run_id <- as.integer(m1[2]); replicate <- 1L
      } else if (length(mN) == 3) {
        run_id <- as.integer(mN[2]); replicate <- as.integer(mN[3])
      } else next

      # FIX: column is base_run_id, not run_id
      if (!(run_id %in% nc_grid$base_run_id)) next

      grid_row <- nc_grid[nc_grid$base_run_id == run_id & nc_grid$simulator == sim &
                           nc_grid$replicate == replicate, ]
      if (nrow(grid_row) != 1) {
        warning(sprintf("Expected exactly 1 null_control_grid match for sim=%s run_id=%d replicate=%d, found %d (%s)",
                         sim, run_id, replicate, nrow(grid_row), f))
        next
      }

      rows[[length(rows) + 1]] <- data.frame(
        file_path = f, method = method_dir, source = sim,
        run_id = run_id, replicate = replicate,
        sparsity = grid_row$sparsity, seed = grid_row$seed,
        stringsAsFactors = FALSE
      )
    }
  }
}

nc_manifest <- do.call(rbind, rows)
cat("\nTotal null-control files discovered:", nrow(nc_manifest), "(expected 270 = 45 x 6)\n")
cat("\nBy method:\n"); print(table(nc_manifest$method))
cat("\nBy source:\n"); print(table(nc_manifest$source))
cat("\nBy replicate:\n"); print(table(nc_manifest$replicate))

write.csv(nc_manifest, "data/processed/step6_3_null_control_manifest.csv", row.names = FALSE)
cat("\nWritten to data/processed/step6_3_null_control_manifest.csv\n")
