# Step 6.4 — Detect Variance Hijacking
#
# Spearman correlation between top 10 PCs and two technical covariates
# (total UMI, batch label). Flag if |rho| > 0.7 for any of top 3 PCs with
# either covariate. Also records the first PC rank correlating more with
# biology (true_group) than with either technical covariate.
#
# Methodological note (flagged before coding, consistent with the Step 6.3
# ARI-vs-random precedent): Spearman correlation assumes an ordinal
# relationship. For simulated data, batch_id is bounded at <=3 levels and
# the grid's own "batch" parameter is framed as an ordinal complexity axis
# (none/simple/complex), so an ordinal treatment is defensible. For real
# data, batch_id is an arbitrary donor/sample identifier (e.g. segerstolpe:
# 10 donor codes; tasic2018: 341 animal IDs) with no natural order, so
# Spearman against an arbitrary integer coding is not statistically
# meaningful there. Both are reported: literal Spearman (per spec, using
# factor-coded batch_id) and the correlation ratio (eta) -- the standard,
# order-independent measure of association between a continuous variable
# and an unordered categorical one (Fisher's correlation ratio; the
# categorical analogue of R^2). "Correlation with biology" (true_group) is
# necessarily measured via eta throughout, since true_group is always
# categorical and never ordinal.
#
# Reuses Step 6.1's UMI cache directly (UMI-per-cell is method-invariant).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args
set.seed(42)

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
work <- manifest[!(manifest$is_null_control %in% TRUE) & manifest$status == "ok", ]
cat(sprintf("Total candidate rows: %d\n", nrow(work)))

umi_cache_dir <- "data/processed/step6_umi_cache"

get_umi <- function(source_file) {
  cache_key <- gsub("[/\\\\]", "_", source_file)
  cache_path <- file.path(umi_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(cache_path)) return(readRDS(cache_path)$umi)
  # Fallback: compute fresh if somehow not cached (shouldn't happen for main-grid files)
  src <- readRDS(source_file)
  Matrix::colSums(assay(src, "counts"))
}

# Correlation ratio (eta): fraction of variance in continuous x explained by
# categorical grouping g. Standard association measure for nominal predictors.
correlation_ratio <- function(x, g) {
  g <- as.character(g)
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NA_real_)
  grand_mean <- mean(x)
  ss_total <- sum((x - grand_mean)^2)
  if (ss_total == 0) return(NA_real_)
  ss_between <- sum(sapply(unique(g), function(lvl) {
    xi <- x[g == lvl]
    length(xi) * (mean(xi) - grand_mean)^2
  }))
  sqrt(ss_between / ss_total)
}

if (dry_run) {
  set.seed(1)
  n_sample <- 500
  large_real <- work[work$source %in% c("pbmc68k", "tabula_sapiens_lung", "tasic2018"), ]
  n_large <- min(9, nrow(large_real))
  large_real_sample <- large_real[sample(nrow(large_real), n_large), ]
  remaining <- work[!(rownames(work) %in% rownames(large_real_sample)), ]
  rest_sample <- remaining[sample(nrow(remaining), n_sample - n_large), ]
  work <- rbind(large_real_sample, rest_sample)
  cat(sprintf("DRY RUN: subset to %d rows (incl. %d large-real rows)\n", nrow(work), n_large))
}

N_PCS <- 10
FLAG_THRESHOLD <- 0.7
TOP_K_FOR_FLAG <- 3

process_row <- function(i) {
  r <- work[i, ]
  out <- list(
    data_type = r$data_type, source = r$source, run_id = r$run_id, method = r$method,
    file_path = r$file_path, n_cells = NA, n_pcs_used = NA,
    max_abs_spearman_umi_top3 = NA, max_abs_spearman_batch_top3 = NA,
    max_eta_batch_top3 = NA, batch_note = NA_character_,
    hijack_flag_literal = NA, hijack_flag_relative = NA,
    first_pc_biology_wins = NA_integer_,
    status = "ok", error_msg = NA_character_
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    tg  <- as.character(obj$true_group)
    bid <- as.character(obj$batch_id)
    n_total <- nrow(emb)
    n_pcs <- min(N_PCS, ncol(emb))

    umi <- get_umi(obj$source_file)
    if (length(umi) != n_total) stop(sprintf("alignment failure: source_file has %d cells, embedding has %d", length(umi), n_total))

    labeled <- !is.na(tg) & tg != ""
    n_batch_levels <- length(unique(bid[labeled]))
    batch_note <- if (n_batch_levels <= 1) "excluded_no_batch_variation" else NA_character_
    bid_int <- if (n_batch_levels > 1) as.integer(factor(bid)) else NULL

    spearman_umi <- numeric(n_pcs)
    spearman_batch <- rep(NA_real_, n_pcs)
    eta_batch <- rep(NA_real_, n_pcs)
    eta_biology <- numeric(n_pcs)

    for (pc in seq_len(n_pcs)) {
      pc_vals <- emb[, pc]
      spearman_umi[pc] <- suppressWarnings(cor(pc_vals, umi, method = "spearman"))
      if (n_batch_levels > 1) {
        spearman_batch[pc] <- suppressWarnings(cor(pc_vals, bid_int, method = "spearman"))
        eta_batch[pc] <- correlation_ratio(pc_vals, bid)
      }
      eta_biology[pc] <- correlation_ratio(pc_vals, tg)
    }

    top3 <- seq_len(min(TOP_K_FOR_FLAG, n_pcs))
    max_umi_top3 <- max(abs(spearman_umi[top3]), na.rm = TRUE)
    max_batch_sp_top3 <- if (n_batch_levels > 1) max(abs(spearman_batch[top3]), na.rm = TRUE) else NA_real_
    max_batch_eta_top3 <- if (n_batch_levels > 1) max(eta_batch[top3], na.rm = TRUE) else NA_real_

    flag_literal <- isTRUE(max_umi_top3 > FLAG_THRESHOLD) ||
                     (!is.na(max_batch_sp_top3) && isTRUE(max_batch_sp_top3 > FLAG_THRESHOLD))
    flag_relative <- isTRUE(max_umi_top3 > FLAG_THRESHOLD) ||
                      (!is.na(max_batch_eta_top3) && isTRUE(max_batch_eta_top3 > FLAG_THRESHOLD))

    # First PC where biology (eta) exceeds BOTH technical measures (UMI |rho|, batch eta)
    technical_strength <- pmax(abs(spearman_umi), ifelse(is.na(eta_batch), 0, eta_batch))
    biology_wins <- which(eta_biology > technical_strength)
    first_win <- if (length(biology_wins) > 0) min(biology_wins) else NA_integer_

    list(n_cells = n_total, n_pcs_used = n_pcs,
         max_abs_spearman_umi_top3 = max_umi_top3,
         max_abs_spearman_batch_top3 = max_batch_sp_top3,
         max_eta_batch_top3 = max_batch_eta_top3,
         batch_note = batch_note,
         hijack_flag_literal = flag_literal, hijack_flag_relative = flag_relative,
         first_pc_biology_wins = first_win,
         status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(n_cells = NA, n_pcs_used = NA, max_abs_spearman_umi_top3 = NA,
         max_abs_spearman_batch_top3 = NA, max_eta_batch_top3 = NA, batch_note = NA_character_,
         hijack_flag_literal = NA, hijack_flag_relative = NA, first_pc_biology_wins = NA_integer_,
         status = "error", error_msg = conditionMessage(e))
  })

  out[names(res)] <- res
  out
}

cat("Starting mclapply (4 workers, mc.preschedule=FALSE)...\n")
t0 <- Sys.time()
results_list <- mclapply(seq_len(nrow(work)), process_row, mc.cores = 4, mc.preschedule = FALSE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

results_df <- do.call(rbind.data.frame, results_list)
n_errors <- sum(results_df$status == "error")
cat(sprintf("Rows processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) print(table(results_df$error_msg))

cat("\nLiteral flag (Spearman only):", sum(results_df$hijack_flag_literal, na.rm=TRUE), "/", nrow(results_df), "\n")
cat("Relative flag (Spearman UMI + eta batch):", sum(results_df$hijack_flag_relative, na.rm=TRUE), "/", nrow(results_df), "\n")
cat("\nFirst-PC-biology-wins distribution:\n")
print(table(results_df$first_pc_biology_wins, useNA = "always"))

suffix <- if (dry_run) "_DRYRUN" else ""
write.csv(results_df, paste0("data/processed/step6_4_variance_hijacking", suffix, ".csv"), row.names = FALSE)
cat("Written.\n")
