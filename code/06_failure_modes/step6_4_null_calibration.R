# Item 6: null-calibration for Step 6.4 (Variance Hijacking), |Spearman|>0.7
# (UMI) and eta>0.7 (batch). Both statistics are correlation-type measures
# whose null distribution depends on n (Spearman) or (n, n_categories)
# (eta) -- not on group/embedding structure -- so bucketing is by n_cells
# alone (Spearman) or (n_cells, n_categories) (eta). Critically, the flag
# logic takes max(|.|) across the top 3 PCs, so the null must do the same
# (3 independent draws per replicate, take the max) to be a fair
# comparison -- this itself bakes in a small multiple-comparisons
# correction within the calibration.

set.seed(42)
N_DRAWS <- 200
TOP_K <- 3

correlation_ratio <- function(x, g) {
  g <- as.character(g)
  grand_mean <- mean(x)
  ss_total <- sum((x - grand_mean)^2)
  if (ss_total == 0) return(NA_real_)
  ss_between <- sum(sapply(unique(g), function(lvl) {
    xi <- x[g == lvl]
    length(xi) * (mean(xi) - grand_mean)^2
  }))
  sqrt(ss_between / ss_total)
}

d <- read.csv("data/processed/step6_4_variance_hijacking.csv", stringsAsFactors = FALSE)
manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)

# ---- Spearman-UMI null: bucket by n_cells alone ----
n_buckets <- sort(unique(d$n_cells[!is.na(d$n_cells)]))
cat("Spearman-UMI calibration buckets:", length(n_buckets), "\n")

calibrate_spearman_max_null <- function(n, n_draws, top_k, seed_base) {
  set.seed(seed_base + n)
  replicate(n_draws, {
    umi_rand <- rnorm(n)
    max_abs <- max(sapply(seq_len(top_k), function(i) {
      abs(suppressWarnings(cor(rnorm(n), umi_rand, method = "spearman")))
    }))
    max_abs
  })
}

spearman_calib <- list()
for (n in n_buckets) {
  null_vals <- calibrate_spearman_max_null(n, N_DRAWS, TOP_K, 3000)
  spearman_calib[[as.character(n)]] <- list(
    n_cells = n, null_mean = mean(null_vals),
    null_95pct = quantile(null_vals, 0.95), null_99pct = quantile(null_vals, 0.99),
    null_max = max(null_vals)
  )
}
spearman_calib_df <- do.call(rbind.data.frame, spearman_calib)
cat("\nSpearman-UMI null calibration (max of top-3 |rho|):\n"); print(spearman_calib_df)

# ---- Eta-batch null: bucket by (n_cells, n_categories) ----
sim_meta <- manifest[manifest$data_type == "simulated", c("source","run_id","method","batch")]
d_sim <- merge(d[d$data_type == "simulated", ], sim_meta, by = c("source","run_id","method"), all.x = TRUE)
sim_buckets <- unique(d_sim[d_sim$batch %in% c("simple","complex") & !is.na(d_sim$n_cells),
                             c("n_cells","batch")])
sim_buckets$n_cat <- ifelse(sim_buckets$batch == "simple", 2, 3)

real_meta <- manifest[manifest$data_type == "real", c("source","n_batches")]
real_meta <- unique(real_meta)
d_real <- merge(d[d$data_type == "real", ], real_meta, by = "source", all.x = TRUE)
real_buckets <- unique(d_real[!is.na(d_real$n_cells) & !is.na(d_real$n_batches) & d_real$n_batches >= 2,
                               c("n_cells","n_batches")])
names(real_buckets)[2] <- "n_cat"

cat("\nEta-batch calibration buckets: simulated =", nrow(sim_buckets), "| real =", nrow(real_buckets), "\n")

calibrate_eta_max_null <- function(n, n_cat, n_draws, top_k, seed_base) {
  set.seed(seed_base + n + n_cat * 97)
  replicate(n_draws, {
    g_rand <- sample(rep(1:n_cat, length.out = n))
    max_eta <- max(sapply(seq_len(top_k), function(i) {
      correlation_ratio(rnorm(n), g_rand)
    }))
    max_eta
  })
}

eta_calib <- list()
for (i in seq_len(nrow(sim_buckets))) {
  n <- sim_buckets$n_cells[i]; ncat <- sim_buckets$n_cat[i]
  null_vals <- calibrate_eta_max_null(n, ncat, N_DRAWS, TOP_K, 4000)
  eta_calib[[paste("sim", n, ncat)]] <- list(
    row_type = "simulated", n_cells = n, n_cat = ncat,
    null_mean = mean(null_vals), null_95pct = quantile(null_vals, 0.95),
    null_99pct = quantile(null_vals, 0.99), null_max = max(null_vals)
  )
}
for (i in seq_len(nrow(real_buckets))) {
  n <- real_buckets$n_cells[i]; ncat <- real_buckets$n_cat[i]
  null_vals <- calibrate_eta_max_null(n, ncat, N_DRAWS, TOP_K, 4000)
  eta_calib[[paste("real", n, ncat)]] <- list(
    row_type = "real", n_cells = n, n_cat = ncat,
    null_mean = mean(null_vals), null_95pct = quantile(null_vals, 0.95),
    null_99pct = quantile(null_vals, 0.99), null_max = max(null_vals)
  )
}
eta_calib_df <- do.call(rbind.data.frame, eta_calib)
cat("\nEta-batch null calibration (max of top-3 eta):\n"); print(eta_calib_df)

# ---- Join back onto main results ----
umi_lookup <- setNames(spearman_calib_df$null_95pct, as.character(spearman_calib_df$n_cells))
d$spearman_umi_null_95pct <- umi_lookup[as.character(d$n_cells)]
d$hijack_flag_umi_calibrated <- d$max_abs_spearman_umi_top3 > d$spearman_umi_null_95pct

sim_eta_lookup <- setNames(
  eta_calib_df$null_95pct[eta_calib_df$row_type == "simulated"],
  paste(eta_calib_df$n_cells[eta_calib_df$row_type == "simulated"], eta_calib_df$n_cat[eta_calib_df$row_type == "simulated"])
)
real_eta_lookup <- setNames(
  eta_calib_df$null_95pct[eta_calib_df$row_type == "real"],
  paste(eta_calib_df$n_cells[eta_calib_df$row_type == "real"], eta_calib_df$n_cat[eta_calib_df$row_type == "real"])
)
d_batch_cat <- merge(d[, c("data_type","source","run_id","method")], sim_meta, by = c("source","run_id","method"), all.x = TRUE)
d_batch_cat$n_cat_sim <- ifelse(d_batch_cat$batch == "simple", 2, ifelse(d_batch_cat$batch == "complex", 3, NA))
d_batch_cat <- merge(d_batch_cat, real_meta, by = "source", all.x = TRUE)
d$n_cat_key <- ifelse(d$data_type == "real", paste(d$n_cells, d_batch_cat$n_batches),
                       paste(d$n_cells, d_batch_cat$n_cat_sim))
d$eta_batch_null_95pct <- ifelse(d$data_type == "real", real_eta_lookup[d$n_cat_key], sim_eta_lookup[d$n_cat_key])
d$hijack_flag_batch_calibrated <- ifelse(is.na(d$eta_batch_null_95pct), NA,
                                          d$max_eta_batch_top3 > d$eta_batch_null_95pct)

d$hijack_flag_calibrated <- mapply(function(u, b) {
  if (isTRUE(u) || isTRUE(b)) TRUE else if (is.na(u) && is.na(b)) NA else FALSE
}, d$hijack_flag_umi_calibrated, d$hijack_flag_batch_calibrated)
d$n_cat_key <- NULL

write.csv(d, "data/processed/step6_4_variance_hijacking_calibrated.csv", row.names = FALSE)
write.csv(spearman_calib_df, "data/processed/step6_4_spearman_umi_null_calibration.csv", row.names = FALSE)
write.csv(eta_calib_df, "data/processed/step6_4_eta_batch_null_calibration.csv", row.names = FALSE)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat("Literal flag (Spearman only, >0.7):", sum(d$hijack_flag_literal, na.rm=TRUE), "\n")
cat("Relative-methodology flag (Spearman UMI + eta batch, >0.7):", sum(d$hijack_flag_relative, na.rm=TRUE), "\n")
cat("Calibrated flag (>95th pct null):", sum(d$hijack_flag_calibrated, na.rm=TRUE), "\n")
