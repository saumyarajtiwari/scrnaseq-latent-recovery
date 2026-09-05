# Item 6: null-calibration for Step 6.1 (Technical Separation), following
# Step 6.3's established precedent -- literal spec flag (AMI>0.5) is kept
# untouched; this adds calibrated context via label-permutation nulls,
# stratified by (n_cells, k_used) for UMI-quartile AMI and additionally by
# batch level (simple/complex) for batch AMI, since n_batch_levels is
# deterministic from the batch grid parameter, not a free variable.

suppressPackageStartupMessages({
  library(aricode)
  library(parallel)
})

set.seed(42)
N_DRAWS <- 200

d <- read.csv("data/processed/step6_1_technical_separation.csv", stringsAsFactors = FALSE)
manifest <- read.csv("data/processed/step3_10_pre_internal_method_backup/embedding_manifest.csv", stringsAsFactors = FALSE)

sim_meta <- manifest[manifest$data_type == "simulated", c("source","run_id","method","batch")]
d <- merge(d, sim_meta, by = c("source","run_id","method"), all.x = TRUE)

# ---- UMI-quartile AMI calibration: bucket by (n_cells, k_used) ----
umi_buckets <- unique(d[!is.na(d$n_cells) & !is.na(d$k_used), c("n_cells","k_used")])
cat("UMI-AMI calibration buckets:", nrow(umi_buckets), "\n")

calibrate_ami_null <- function(n, k1, k2, n_draws) {
  set.seed(1000 + n + k1 * 7 + k2 * 13)
  replicate(n_draws, {
    a <- sample(rep(1:k1, length.out = n))
    b <- sample(rep(1:k2, length.out = n))
    aricode::AMI(a, b)
  })
}

umi_calib <- list()
for (i in seq_len(nrow(umi_buckets))) {
  n <- umi_buckets$n_cells[i]; k <- umi_buckets$k_used[i]
  ami_null <- calibrate_ami_null(n, k, 4, N_DRAWS)  # UMI quartile = 4 levels
  umi_calib[[paste(n, k)]] <- list(
    n_cells = n, k_used = k,
    null_mean = mean(ami_null), null_95pct = quantile(ami_null, 0.95),
    null_99pct = quantile(ami_null, 0.99), null_max = max(ami_null)
  )
}
umi_calib_df <- do.call(rbind.data.frame, umi_calib)
cat("\nUMI-AMI null calibration summary:\n"); print(umi_calib_df)

# ---- Batch AMI calibration: bucket by (n_cells, k_used, n_batch_levels) ----
# SIMULATED rows: n_batch_levels inferred empirically from a representative
# file per (n_cells,k,batch grid-param) bucket (simple/complex only; none
# is excluded upstream, has no batch variation).
# REAL rows: n_batch_levels comes directly from d$n_batches, already
# populated for exactly these 36 rows from real_data_inventory.csv (Step
# 6.1's original script computes real batch AMI from actual batch_id, not
# a grid parameter -- using the same real cardinality here for parity).
batch_rows_sim <- d[d$batch %in% c("simple","complex") & !is.na(d$n_cells) & !is.na(d$k_used), ]
batch_buckets_sim <- unique(batch_rows_sim[, c("n_cells","k_used","batch")])
cat("\nBatch-AMI calibration buckets (simulated):", nrow(batch_buckets_sim), "\n")

get_n_batch_levels <- function(n_cells, k_used, batch_val) {
  rep_row <- batch_rows_sim[batch_rows_sim$n_cells == n_cells & batch_rows_sim$k_used == k_used &
                             batch_rows_sim$batch == batch_val, ][1, ]
  obj <- tryCatch(readRDS(rep_row$file_path), error = function(e) NULL)
  if (is.null(obj)) return(NA_integer_)
  length(unique(obj$batch_id[!is.na(obj$true_group) & obj$true_group != ""]))
}

batch_calib <- list()
for (i in seq_len(nrow(batch_buckets_sim))) {
  n <- batch_buckets_sim$n_cells[i]; k <- batch_buckets_sim$k_used[i]; bval <- batch_buckets_sim$batch[i]
  n_batch_lvl <- get_n_batch_levels(n, k, bval)
  if (is.na(n_batch_lvl) || n_batch_lvl < 2) next
  ami_null <- calibrate_ami_null(n, k, n_batch_lvl, N_DRAWS)
  batch_calib[[paste("sim", n, k, bval)]] <- list(
    row_type = "simulated", n_cells = n, k_used = k, batch = bval, n_batch_levels = n_batch_lvl,
    null_mean = mean(ami_null), null_95pct = quantile(ami_null, 0.95),
    null_99pct = quantile(ami_null, 0.99), null_max = max(ami_null)
  )
}

# Real-data buckets: use d$n_batches directly (already correct per-row).
real_rows <- d[d$data_type == "real" & !is.na(d$n_cells) & !is.na(d$k_used) &
               !is.na(d$n_batches) & d$n_batches >= 2, ]
real_buckets <- unique(real_rows[, c("n_cells","k_used","n_batches")])
cat("Batch-AMI calibration buckets (real):", nrow(real_buckets), "\n")
for (i in seq_len(nrow(real_buckets))) {
  n <- real_buckets$n_cells[i]; k <- real_buckets$k_used[i]; nb <- real_buckets$n_batches[i]
  ami_null <- calibrate_ami_null(n, k, nb, N_DRAWS)
  batch_calib[[paste("real", n, k, nb)]] <- list(
    row_type = "real", n_cells = n, k_used = k, batch = NA_character_, n_batch_levels = nb,
    null_mean = mean(ami_null), null_95pct = quantile(ami_null, 0.95),
    null_99pct = quantile(ami_null, 0.99), null_max = max(ami_null)
  )
}
batch_calib_df <- do.call(rbind.data.frame, batch_calib)
cat("\nBatch-AMI null calibration summary:\n"); print(batch_calib_df)

# ---- Join calibration back onto main results, add relative-flag columns ----
d$umi_key <- paste(d$n_cells, d$k_used)
umi_lookup <- setNames(umi_calib_df$null_95pct, paste(umi_calib_df$n_cells, umi_calib_df$k_used))
d$umi_ami_null_95pct <- umi_lookup[d$umi_key]
d$technical_separation_flag_umi_relative <- d$umi_quartile_ami > d$umi_ami_null_95pct

# Simulated rows key on (n_cells, k_used, batch grid-param); real rows key
# on (n_cells, k_used, n_batches) since they have no `batch` grid parameter.
sim_batch_lookup <- setNames(
  batch_calib_df$null_95pct[batch_calib_df$row_type == "simulated"],
  paste(batch_calib_df$n_cells[batch_calib_df$row_type == "simulated"],
        batch_calib_df$k_used[batch_calib_df$row_type == "simulated"],
        batch_calib_df$batch[batch_calib_df$row_type == "simulated"])
)
real_batch_lookup <- setNames(
  batch_calib_df$null_95pct[batch_calib_df$row_type == "real"],
  paste(batch_calib_df$n_cells[batch_calib_df$row_type == "real"],
        batch_calib_df$k_used[batch_calib_df$row_type == "real"],
        batch_calib_df$n_batch_levels[batch_calib_df$row_type == "real"])
)
d$batch_ami_null_95pct <- ifelse(
  d$data_type == "real",
  real_batch_lookup[paste(d$n_cells, d$k_used, d$n_batches)],
  sim_batch_lookup[paste(d$n_cells, d$k_used, d$batch)]
)
d$technical_separation_flag_batch_relative <- ifelse(
  is.na(d$batch_ami_null_95pct), NA,
  d$batch_ami > d$batch_ami_null_95pct
)

d$technical_separation_flag_relative <- mapply(function(u, b) {
  if (isTRUE(u) || isTRUE(b)) TRUE else if (is.na(u) && is.na(b)) NA else FALSE
}, d$umi_quartile_ami > d$umi_ami_null_95pct, d$technical_separation_flag_batch_relative)

d$umi_key <- NULL; d$batch_key <- NULL

write.csv(d, "data/processed/step6_1_technical_separation_calibrated.csv", row.names = FALSE)
write.csv(umi_calib_df, "data/processed/step6_1_umi_ami_null_calibration.csv", row.names = FALSE)
write.csv(batch_calib_df, "data/processed/step6_1_batch_ami_null_calibration.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Literal flag (AMI>0.5) count:", sum(d$technical_separation_flag, na.rm=TRUE), "\n")
cat("Relative flag (>95th pct null) count:", sum(d$technical_separation_flag_relative, na.rm=TRUE), "\n")
