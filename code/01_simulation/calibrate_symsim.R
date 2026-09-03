# =============================================================================
# calibrate_symsim.R
# Phase A calibration for SymSim simulation runs.
#
# Tasks:
#   A3: Sigma -> actual sparsity (at baseline depth and dropout)
#   A4: depth_mean -> actual mean UMI/cell (per alpha_mean level)
#   A5: alpha_mean -> dropout separation confirmation
#
# Key API facts confirmed:
#   - true_counts <- true_result[["counts"]]
#   - cell_meta   <- true_result[["cell_meta"]]  (pop column = cell type label)
#   - True2ObservedCounts has no seed arg -> set.seed(run_id) before each call
#   - depth_mean does NOT map linearly to UMI/cell (empirical calibration needed)
#   - alpha_mean = capture efficiency (higher = more UMI, less dropout)
#   - bimod fixed at 0 throughout (no expression bimodality added)
#   - ngenes = 2000 to match actual simulation
#
# Output:
#   data/simulated/symsim_calib_sigma.csv
#   data/simulated/symsim_calib_depth.csv
# =============================================================================

suppressPackageStartupMessages(library(SymSim))

N_GENES   <- 2000L
N_CELLS   <- 200L
N_REPS    <- 2L

cat("=== SymSim Calibration ===\n")
cat(sprintf("ngenes=%d, ncells=%d, reps=%d\n", N_GENES, N_CELLS, N_REPS))
cat(sprintf("Started: %s\n\n", format(Sys.time())))

# -----------------------------------------------------------------------------
# Helper: one SimulateTrueCounts + True2ObservedCounts call
# seed used for both randseed and set.seed before True2ObservedCounts
# -----------------------------------------------------------------------------
run_one <- function(seed, sigma, depth_mean, alpha_mean) {
  tr <- SimulateTrueCounts(
    ncells_total = N_CELLS,
    min_popsize  = N_CELLS,
    i_minpop     = 1L,
    ngenes       = N_GENES,
    evf_type     = "one.population",
    randseed     = seed,
    Sigma        = sigma,
    n_de_evf     = 8L,
    bimod        = 0,
    vary         = "s"
  )
  true_counts <- tr[["counts"]]
  cell_meta   <- tr[["cell_meta"]]

  set.seed(seed)
  obs <- True2ObservedCounts(
    true_counts = true_counts,
    meta_cell   = cell_meta,
    protocol    = "UMI",
    alpha_mean  = alpha_mean,
    alpha_sd    = 0.01,
    gene_len    = rep(1000L, N_GENES),
    depth_mean  = depth_mean,
    depth_sd    = max(50L, round(depth_mean * 0.10))
  )
  obs_counts <- obs[[1]]

  list(
    mean_umi = round(mean(colSums(obs_counts)), 2),
    sparsity = round(sum(obs_counts == 0L) / length(obs_counts), 4)
  )
}

# =============================================================================
# A3 — Sigma calibration: Sigma -> sparsity
# Fixed: alpha_mean=0.08 (low dropout), depth_mean=5000 (moderate)
# =============================================================================
cat("PART A3 — Sigma -> sparsity\n")
cat("    [alpha_mean=0.08, depth_mean=5000]\n\n")

sigma_vals <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)
sigma_df   <- data.frame(sigma=sigma_vals,
                          mean_sparsity=NA_real_,
                          mean_umi=NA_real_)

for (i in seq_along(sigma_vals)) {
  s    <- sigma_vals[i]
  reps <- sapply(seq_len(N_REPS), function(r) {
    res <- tryCatch(run_one(r * 100L, s, 5000, 0.08),
                    error = function(e) list(sparsity=NA, mean_umi=NA))
    c(res$sparsity, res$mean_umi)
  })
  sigma_df$mean_sparsity[i] <- round(mean(reps[1,], na.rm=TRUE), 4)
  sigma_df$mean_umi[i]      <- round(mean(reps[2,], na.rm=TRUE), 1)
  cat(sprintf("  Sigma=%.1f -> sparsity=%.4f  UMI=%.1f\n",
              s, sigma_df$mean_sparsity[i], sigma_df$mean_umi[i]))
}

cat("\n  Target -> recommended Sigma:\n")
for (t in c(0.90, 0.93, 0.95, 0.97, 0.99)) {
  idx <- which.min(abs(sigma_df$mean_sparsity - t))
  cat(sprintf("  target %.2f -> Sigma=%.1f (actual=%.4f)\n",
              t, sigma_df$sigma[idx], sigma_df$mean_sparsity[idx]))
}

# =============================================================================
# A4 — Depth calibration: depth_mean -> actual UMI/cell per alpha_mean
# =============================================================================
cat("\nPART A4 — depth_mean -> actual UMI/cell\n")
cat("    [Sigma=0.4, per alpha_mean level]\n\n")

depth_vals <- c(200, 500, 1000, 2000, 5000, 10000, 25000, 50000)
alpha_levels <- c(
  none = 0.97,
  low  = 0.08,
  high = 0.02
)

depth_rows <- list()
for (aname in names(alpha_levels)) {
  alpha <- alpha_levels[aname]
  cat(sprintf("  alpha_mean=%.2f (%s dropout):\n", alpha, aname))
  for (dm in depth_vals) {
    res <- tryCatch(run_one(42L, 0.4, dm, alpha),
                    error = function(e) list(mean_umi=NA, sparsity=NA))
    cat(sprintf("    depth_mean=%6d -> UMI=%.1f  sparsity=%.4f\n",
                dm, res$mean_umi, res$sparsity))
    depth_rows[[length(depth_rows)+1]] <- data.frame(
      dropout_label = aname,
      alpha_mean    = alpha,
      depth_mean    = dm,
      actual_umi    = res$mean_umi,
      actual_sparsity = res$sparsity
    )
  }
  cat(sprintf("\n    -> recommended depth_mean for targets:\n"))
  umi_results <- sapply(depth_vals, function(dm) {
    tryCatch(run_one(43L, 0.4, dm, alpha)$mean_umi, error=function(e) NA)
  })
  for (target in c(500, 2000, 10000)) {
    idx <- which.min(abs(umi_results - target))
    cat(sprintf("      target %5d UMI -> depth_mean=%6d (actual=%.1f)\n",
                target, depth_vals[idx], umi_results[idx]))
  }
  cat("\n")
}

depth_df <- do.call(rbind, depth_rows)

# =============================================================================
# A5 — Dropout separation confirmation
# =============================================================================
cat("PART A5 — alpha_mean dropout separation\n")
cat("    [Sigma=0.4, depth_mean=5000]\n\n")

for (aname in names(alpha_levels)) {
  alpha <- alpha_levels[aname]
  res   <- tryCatch(run_one(42L, 0.4, 5000, alpha),
                    error=function(e) list(mean_umi=NA, sparsity=NA))
  cat(sprintf("  %-4s (alpha=%.2f): UMI=%.1f  sparsity=%.4f\n",
              aname, alpha, res$mean_umi, res$sparsity))
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================
cat("\nWriting calibration tables...\n")
write.csv(sigma_df, "data/simulated/symsim_calib_sigma.csv", row.names=FALSE)
write.csv(depth_df, "data/simulated/symsim_calib_depth.csv", row.names=FALSE)
cat("  data/simulated/symsim_calib_sigma.csv\n")
cat("  data/simulated/symsim_calib_depth.csv\n")
cat(sprintf("\n=== Calibration complete: %s ===\n", format(Sys.time())))
cat("Review output above before updating param_dict.R\n")
