# Step 6.8 — Cross-Reference Failure Modes with Phase Boundaries
#
# Maps each failure-mode's flag rate (Steps 6.1-6.4, 6.6-6.7) onto Step 5.8's
# geometric failure-boundary thresholds (interpolated stressor value at
# which subspace_recovery_score crosses 0.5), per (axis, simulator, method).
# Confirms whether independently-derived failure-mode evidence transitions
# at the same stressor regime as the geometric recovery-score boundary.
#
# Scope notes:
# - 24/108 Step 5.8 rows are already flagged unreliable (SymSim
#   batch/dropout anomalies, clipping inertness for scdesign3/SymSim) --
#   excluded from the core alignment analysis (comparing against an
#   already-known-unreliable threshold is not meaningful either way),
#   reported separately.
# - Step 6.3 (Phantom Clustering) only applies to null-controls, which vary
#   ONLY sparsity -- included for the sparsity axis only.
# - Step 6.5 (Over-Smoothing) excluded: the pancreas DPT branch is not
#   grid-axis-based (real data, no simulated stressor axes), and the
#   simulated variance-ratio branch's scope (log2FC bands, 4 methods,
#   0/31,640 flagged) is too narrow/orthogonal to map cleanly onto these 6
#   axes.
# - Marginal aggregation: flag rate per axis level is averaged across all
#   other grid dimensions (this is a full-factorial design, not one-axis-
#   at-a-time). NOTE: this does NOT match Step 5's axis-sweep methodology --
#   step5_1_organize_by_axis.R holds every other axis at a single fixed
#   baseline (OFAT design), so Step 5.8's per-axis thresholds characterize
#   behavior at ONE grid point, while this script's flag rates characterize
#   the FULL marginal grid (e.g. 3,643 rows vs. 1, for splatter/shiftedlog/
#   batch=complex). This scope mismatch is the resolved explanation for the
#   3 combinations flagged as an open contradiction in
#   docs/step6_9_failure_mode_review.md -- see that document and
#   PROJECT_HANDOVER.md for the full writeup.

AXIS_LEVEL_ORDER <- list(
  sparsity = c("0.7","0.8","0.9","0.95","0.98"),
  depth = c("500","2000","10000"),
  dropout = c("none","low","high"),
  separability = c("high","medium","low"),
  batch = c("none","simple","complex"),
  clipping = c("none","clip99","log_stabilized")
)

thresholds <- read.csv("data/processed/step5_8_critical_thresholds.csv", stringsAsFactors = FALSE)
reliable_thresholds <- thresholds[thresholds$reliable == TRUE, ]
unreliable_thresholds <- thresholds[thresholds$reliable == FALSE, ]
cat(sprintf("Reliable threshold rows: %d | Unreliable (excluded from core analysis): %d\n",
            nrow(reliable_thresholds), nrow(unreliable_thresholds)))

cat("Loading failure-mode results...\n")
d61 <- read.csv("data/processed/step6_1_technical_separation.csv", stringsAsFactors = FALSE)
d62 <- read.csv("data/processed/step6_2_cluster_collapse.csv", stringsAsFactors = FALSE)
d63 <- read.csv("data/processed/step6_3_phantom_clustering.csv", stringsAsFactors = FALSE)
d64 <- read.csv("data/processed/step6_4_variance_hijacking.csv", stringsAsFactors = FALSE)
d66 <- read.csv("data/processed/step6_6_neighborhood_collapse.csv", stringsAsFactors = FALSE)
d67 <- read.csv("data/processed/step6_7_subspace_rotation_slippage.csv", stringsAsFactors = FALSE)

# Attach grid-axis columns to failure-mode tables that don't already carry
# them (6.1/6.2/6.4/6.6 only have source/run_id/method -- join against the
# manifest for axis values). 6.3 and 6.7 already carry axis columns directly.
manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
axis_lookup <- manifest[manifest$data_type == "simulated" & !(manifest$is_null_control %in% TRUE),
                         c("source","run_id","method","sparsity","depth","dropout",
                           "separability","batch","clipping")]
axis_lookup <- unique(axis_lookup)
join_axes <- function(df) {
  df <- df[df$data_type == "simulated", ]
  merge(df, axis_lookup, by = c("source","run_id","method"), all.x = TRUE)
}

d61s <- join_axes(d61); d61s$flag_col <- d61s$technical_separation_flag
d62s <- join_axes(d62); d62s$flag_col <- d62s$n_pairs_collapsed > 0
d64s <- join_axes(d64); d64s$flag_col <- d64s$hijack_flag_relative
d66s <- join_axes(d66); d66s$flag_col <- d66s$overlap_flag
d67$flag_col <- d67$flag_pc1  # already has axis columns; using PC1 (dominant slip axis)
d63$flag_col <- d63$phantom_flag_ch_relative  # already has sparsity column

failure_tables <- list(
  step6_1_technical_separation = d61s,
  step6_2_cluster_collapse = d62s,
  step6_3_phantom_clustering = d63,
  step6_4_variance_hijacking = d64s,
  step6_6_neighborhood_collapse = d66s,
  step6_7_subspace_rotation_slippage = d67
)

get_marginal_flag_rate <- function(df, axis_col, simulator, method) {
  sub <- df[df$source == simulator & df$method == method, ]
  if (!(axis_col %in% colnames(sub)) || nrow(sub) == 0) return(NULL)
  sub <- sub[!is.na(sub[[axis_col]]) & !is.na(sub$flag_col), ]
  if (nrow(sub) == 0) return(NULL)
  agg <- aggregate(flag_col ~ get(axis_col), data = sub, FUN = mean)
  colnames(agg) <- c("level", "flag_rate")
  agg
}

results <- list()

for (i in seq_len(nrow(reliable_thresholds))) {
  r <- reliable_thresholds[i, ]
  axis_col <- r$axis
  level_order <- AXIS_LEVEL_ORDER[[axis_col]]

  applicable_tables <- names(failure_tables)
  if (axis_col != "sparsity") applicable_tables <- setdiff(applicable_tables, "step6_3_phantom_clustering")

  for (tbl_name in applicable_tables) {
    df <- failure_tables[[tbl_name]]
    agg <- get_marginal_flag_rate(df, axis_col, r$simulator, r$method)
    if (is.null(agg) || nrow(agg) < 2) next

    agg$level <- factor(as.character(agg$level), levels = level_order, ordered = TRUE)
    agg <- agg[!is.na(agg$level), ]
    agg <- agg[order(agg$level), ]
    if (nrow(agg) < 2) next

    # Alignment check
    alignment <- NA_character_
    if (r$threshold_status == "always_below_threshold") {
      # Recovery always fails -> expect uniformly HIGH flag rate
      alignment <- if (mean(agg$flag_rate) > 0.5) "consistent_uniform_high" else "inconsistent_low_flag_despite_always_failing"
    } else if (r$threshold_status == "always_above_threshold") {
      # Recovery never fails -> expect uniformly LOW flag rate
      alignment <- if (mean(agg$flag_rate) < 0.5) "consistent_uniform_low" else "inconsistent_high_flag_despite_never_failing"
    } else {
      # interpolated or between_X_and_Y: expect flag rate to be lower before
      # threshold_value's corresponding level and higher after
      levels_numeric <- as.numeric(agg$level)
      below_half <- levels_numeric <= median(levels_numeric)
      mean_before <- mean(agg$flag_rate[below_half])
      mean_after <- mean(agg$flag_rate[!below_half])
      if (is.nan(mean_before) || is.nan(mean_after)) {
        alignment <- "insufficient_levels"
      } else if (mean_after > mean_before) {
        alignment <- "consistent_increasing"
      } else {
        alignment <- "inconsistent_not_increasing"
      }
    }

    results[[length(results) + 1]] <- data.frame(
      axis = axis_col, simulator = r$simulator, method = r$method,
      method_family = r$method_family,
      failure_mode = tbl_name, threshold_status = r$threshold_status,
      threshold_value = r$threshold_value, mean_flag_rate = mean(agg$flag_rate),
      min_flag_rate = min(agg$flag_rate), max_flag_rate = max(agg$flag_rate),
      alignment = alignment,
      ground_truth_caveat = r$ground_truth_caveat
    )
  }
}

results_df <- do.call(rbind, results)
write.csv(results_df, "data/processed/step6_8_cross_reference.csv", row.names = FALSE)
write.csv(unreliable_thresholds, "data/processed/step6_8_excluded_unreliable_thresholds.csv", row.names = FALSE)

cat("\nTotal comparison rows:", nrow(results_df), "\n")
cat("\nAlignment summary:\n")
print(table(results_df$alignment))
cat("\nAlignment by failure mode:\n")
print(table(results_df$failure_mode, results_df$alignment))
cat("\nWritten to data/processed/step6_8_cross_reference.csv\n")
