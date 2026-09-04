d <- read.csv("data/processed/step5_axis_sweep_table.csv", stringsAsFactors = FALSE)

METHOD_LABELS <- c(
  pca_raw = "Raw PCA", pca_libnorm = "Libnorm PCA", pca_log = "Log-PCA",
  pca_shiftedlog = "Shifted Log-PCA", pca_sctransform_v2 = "Pearson Residual PCA (SCT v2)",
  pca_glmpca = "GLM-PCA"
)

# Item 1 fix: ground-truth-linearity stratification (see PROJECT_HANDOVER.md
# ground-truth definition section, and param_dict.R's method_family lookup).
# subspace_recovery_score is computed against a LINEAR ground truth for every
# method; for the two nonlinear methods this partly measures "how linear is
# this method's recovered subspace," not purely biological recovery. No rows
# are removed -- nonlinear-method rows just carry an explicit caveat flag so
# they are never silently pooled with the four linear methods downstream.
METHOD_FAMILY <- c(
  pca_raw = "linear", pca_libnorm = "linear", pca_log = "linear",
  pca_shiftedlog = "linear", pca_sctransform_v2 = "nonlinear", pca_glmpca = "nonlinear"
)
GT_LINEARITY_CAVEAT <- paste(
  "subspace_recovery_score computed against a linearly-defined ground truth;",
  "for this nonlinear method the score partly reflects representational",
  "mismatch, not purely biological recovery -- see PROJECT_HANDOVER.md"
)

# Reliability findings from Step 5.2-5.7's investigation, baked in explicitly
RELIABILITY <- list(
  sparsity     = c(scdesign3=TRUE,  splatter=TRUE, symsim=TRUE),
  depth        = c(scdesign3=TRUE,  splatter=TRUE, symsim=TRUE),
  separability = c(scdesign3=TRUE,  splatter=TRUE, symsim=TRUE),
  dropout      = c(scdesign3=TRUE,  splatter=TRUE, symsim=FALSE), # symsim non-monotonic, not swap-fixable
  clipping     = c(scdesign3=FALSE, splatter=TRUE, symsim=FALSE), # confirmed inert for both
  batch        = c(scdesign3=TRUE,  splatter=TRUE, symsim=FALSE)  # symsim shows unexplained anomaly
)
RELIABILITY_NOTE <- list(
  sparsity=NA, depth=NA, separability=NA,
  dropout="symsim excluded: dropout non-monotonic, confirmed not a simple label swap (see docs/step5_axis_findings.md)",
  clipping="scdesign3 and symsim excluded: clipping confirmed generatively inert for both (see docs/step5_axis_findings.md)",
  batch="symsim result reported but flagged unreliable: unexplained anomalous positive relationship, confounding and inertness ruled out (see docs/step5_axis_findings.md)"
)

NUMERIC_AXES <- list(sparsity = "linear", depth = "log10")
CATEGORICAL_ORDER <- list(
  separability = c("low","medium","high"),
  batch        = c("none","simple","complex"),
  dropout      = c("none","low","high"),
  clipping     = c("none","log_stabilized","clip99")
)

interpolate_threshold <- function(x, y, scale = "linear") {
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  if (scale == "log10") x <- log10(x)
  if (all(y >= 0.5)) return(list(value = NA, status = "always_above_threshold"))
  if (all(y < 0.5)) return(list(value = NA, status = "always_below_threshold"))
  for (i in seq_len(length(x) - 1)) {
    if ((y[i] >= 0.5 && y[i+1] < 0.5) || (y[i] < 0.5 && y[i+1] >= 0.5)) {
      frac <- (0.5 - y[i]) / (y[i+1] - y[i])
      x_thresh <- x[i] + frac * (x[i+1] - x[i])
      if (scale == "log10") x_thresh <- 10^x_thresh
      return(list(value = x_thresh, status = "interpolated"))
    }
  }
  list(value = NA, status = "no_single_crossing")
}

categorical_transition <- function(levels_order, y) {
  y <- y[match(levels_order, names(y))]
  if (all(y >= 0.5)) return("always_above_threshold")
  if (all(y < 0.5)) return("always_below_threshold")
  for (i in seq_len(length(levels_order) - 1)) {
    if ((y[i] >= 0.5 && y[i+1] < 0.5) || (y[i] < 0.5 && y[i+1] >= 0.5)) {
      return(paste0("between_", levels_order[i], "_and_", levels_order[i+1]))
    }
  }
  "no_single_crossing"
}

results <- list()

for (axis in names(NUMERIC_AXES)) {
  flag_col <- paste0("in_", axis, "_sweep")
  sub <- d[d[[flag_col]] == TRUE, ]
  for (sim in unique(sub$source)) {
    for (method in unique(sub$method)) {
      sel <- sub[sub$source == sim & sub$method == method, ]
      r <- interpolate_threshold(sel[[axis]], sel$subspace_recovery_score, scale = NUMERIC_AXES[[axis]])
      results[[length(results)+1]] <- data.frame(
        axis = axis, simulator = sim, method = method, method_label = METHOD_LABELS[method],
        method_family = METHOD_FAMILY[method],
        threshold_value = r$value, threshold_status = r$status,
        reliable = RELIABILITY[[axis]][sim], reliability_note = RELIABILITY_NOTE[[axis]],
        ground_truth_caveat = if (METHOD_FAMILY[method] == "nonlinear") GT_LINEARITY_CAVEAT else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
}

for (axis in names(CATEGORICAL_ORDER)) {
  flag_col <- paste0("in_", axis, "_sweep")
  sub <- d[d[[flag_col]] == TRUE, ]
  for (sim in unique(sub$source)) {
    for (method in unique(sub$method)) {
      sel <- sub[sub$source == sim & sub$method == method, ]
      y <- setNames(sel$subspace_recovery_score, sel[[axis]])
      status <- categorical_transition(CATEGORICAL_ORDER[[axis]], y)
      results[[length(results)+1]] <- data.frame(
        axis = axis, simulator = sim, method = method, method_label = METHOD_LABELS[method],
        method_family = METHOD_FAMILY[method],
        threshold_value = NA, threshold_status = status,
        reliable = RELIABILITY[[axis]][sim], reliability_note = RELIABILITY_NOTE[[axis]],
        ground_truth_caveat = if (METHOD_FAMILY[method] == "nonlinear") GT_LINEARITY_CAVEAT else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
}

final <- do.call(rbind, results)
OUT_PATH <- "data/processed/step5_8_critical_thresholds.csv"
write.csv(final, OUT_PATH, row.names = FALSE)

cat("=== STEP 5.8 COMPLETE ===\n")
cat("Total threshold rows:", nrow(final), "\n\n")
cat("Reliable rows only, by axis:\n")
print(table(final$axis[final$reliable]))
cat("\nThreshold status distribution:\n")
print(table(final$threshold_status))
cat("\nRows carrying ground-truth-linearity caveat (nonlinear methods):",
    sum(!is.na(final$ground_truth_caveat)), "of", nrow(final), "\n")
