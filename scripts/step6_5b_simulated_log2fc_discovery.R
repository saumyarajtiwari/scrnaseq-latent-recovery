# Step 6.5b — discover simulated rows near log2FC = 0.25 / 0.5 bands, using
# the same formula as the original EDA 1.3 investigation (log2(mean+1)
# pseudocount, all group-pairs per file, mean_abs_log2fc summary). Uses the
# already-stored true_group_means metadata directly (no per-cell recomputation
# needed), so this is much cheaper than the original 32,814-file EDA pass.

suppressPackageStartupMessages(library(parallel))

sim_files <- list.files("data/simulated/sce", pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
sim_files <- sim_files[!grepl("_rep[23]\\.rds$", sim_files)]  # exclude null-control replicates (n_groups=1, no log2FC meaningful)
cat("Total simulated files to scan:", length(sim_files), "\n")

process_file <- function(fp) {
  res <- tryCatch({
    obj <- readRDS(fp)
    tgm <- S4Vectors::metadata(obj)$true_group_means
    if (is.null(tgm) || ncol(tgm) < 2) stop("no true_group_means or <2 groups")

    log_tgm <- log2(tgm + 1)
    n_groups <- ncol(tgm)
    pairs <- combn(n_groups, 2, simplify = FALSE)
    abs_fc <- sapply(pairs, function(p) abs(log_tgm[, p[1]] - log_tgm[, p[2]]))

    list(mean_abs_log2fc = mean(abs_fc), p90_abs_log2fc = quantile(abs_fc, 0.9, names = FALSE),
         n_groups = n_groups, status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(mean_abs_log2fc = NA, p90_abs_log2fc = NA, n_groups = NA,
         status = "error", error_msg = conditionMessage(e))
  })
  res$file_path <- fp
  res
}

cat("Processing (4 workers)...\n")
t0 <- Sys.time()
results_list <- mclapply(sim_files, process_file, mc.cores = 4, mc.preschedule = TRUE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

results_df <- do.call(rbind.data.frame, results_list)
n_errors <- sum(results_df$status == "error")
cat(sprintf("Processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) print(table(results_df$error_msg))

cat("\nmean_abs_log2fc distribution:\n")
print(summary(results_df$mean_abs_log2fc))

write.csv(results_df, "data/processed/step6_5b_simulated_log2fc.csv", row.names = FALSE)
cat("\nWritten to data/processed/step6_5b_simulated_log2fc.csv\n")

for (target in c(0.25, 0.5)) {
  band <- results_df[!is.na(results_df$mean_abs_log2fc) &
                      abs(results_df$mean_abs_log2fc - target) <= 0.05, ]
  cat(sprintf("\nRows within +/-0.05 of log2FC=%.2f: %d\n", target, nrow(band)))
}
