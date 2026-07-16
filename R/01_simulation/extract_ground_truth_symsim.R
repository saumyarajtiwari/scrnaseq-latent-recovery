# =============================================================================
# extract_ground_truth_symsim.R
#
# Extracts per-fit_key "true" biological signal for all 244 unique SymSim
# fit_keys: group-mean expression from SimulateTrueCounts() output, BEFORE
# True2ObservedCounts() adds technical noise. Uses identical batch/seed
# construction as production simulate_symsim.R (same seed = same cells).
#
# Every fit_key's population assignment is validated by exact match
# (identical()) against the corresponding already-validated production
# file's cell_meta$true_group -- not a sample check, all 244. This is
# possible because post-hoc masking (applied downstream in production)
# only zeroes count values; it never alters cell-to-population identity.
#
# Genes are synthetic (Gene1..Gene2000); raw tr$counts carries no native
# rownames (confirmed separately), so GeneN naming is applied explicitly
# here rather than assumed.
#
# Read-only w.r.t. param_grid.csv, param_dict.R, and all of
# data/simulated/symsim/*.rds (read only for validation, never modified).
# Writes only to:
#   data/simulated/ground_truth/symsim/symsim_truth_run_NNNNN.rds
#   data/simulated/ground_truth/symsim_manifest.csv
# =============================================================================

suppressPackageStartupMessages({
  library(SymSim)
  library(Matrix)
})

source("R/01_simulation/param_dict.R")
param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors=FALSE)

cat("=== SymSim Ground-Truth Extraction ===\n")

fit_cols <- c("separability","n_cells","batch","dropout","depth")
param_grid$fit_key <- apply(param_grid[, fit_cols], 1, function(r) paste(r, collapse="_"))
unique_fits <- unique(param_grid$fit_key)
cat(sprintf("Unique fit_keys: %d\n\n", length(unique_fits)))

N_GENES <- 2000L
SYMSIM_SIGMA <- 0.4
OUT_DIR <- "data/simulated/ground_truth/symsim"
PROD_DIR <- "data/simulated/symsim"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

manifest <- data.frame(fit_key=character(), representative_run_id=integer(),
                        n_populations=integer(), n_cells_total=integer(),
                        matches_production_pop=logical(), output_path=character(),
                        stringsAsFactors=FALSE)

n_mismatch <- 0L
t_start <- Sys.time()

for (i in seq_along(unique_fits)) {
  fk <- unique_fits[i]
  row1   <- param_grid[param_grid$fit_key == fk, ][1, ]
  run_id <- as.integer(row1$run_id)
  out_path <- file.path(OUT_DIR, sprintf("symsim_truth_run_%05d.rds", run_id))

  if (file.exists(out_path)) {
    cat(sprintf("[%3d/%d] SKIP (exists) fit_key=%s\n", i, length(unique_fits), fk))
    manifest <- rbind(manifest, data.frame(fit_key=fk, representative_run_id=run_id,
                       n_populations=NA, n_cells_total=NA, matches_production_pop=NA,
                       output_path=out_path))
    next
  }

  sep_p    <- symsim_separability[[row1$separability]]
  n_cells  <- as.integer(row1$n_cells)
  is_null  <- isTRUE(as.logical(row1$is_null_control)) || row1$separability == "null"
  evf_type <- sep_p$evf_type
  n_de_evf <- if (is_null) 8L else as.integer(sep_p$n_de_evf)

  if (row1$batch == "none") {
    batch_sizes <- n_cells; batch_seeds <- run_id
  } else if (row1$batch == "simple") {
    b1 <- floor(n_cells/2L)
    batch_sizes <- c(b1, n_cells-b1); batch_seeds <- c(run_id+10000L, run_id+20000L)
  } else {
    b1 <- round(n_cells*0.50); b2 <- round(n_cells*0.30)
    batch_sizes <- c(b1, b2, n_cells-b1-b2)
    batch_seeds <- c(run_id+10000L, run_id+20000L, run_id+30000L)
  }

  all_counts <- vector("list", length(batch_sizes))
  all_pops   <- vector("list", length(batch_sizes))

  for (b in seq_along(batch_sizes)) {
    bs    <- as.integer(batch_sizes[b])
    bseed <- as.integer(batch_seeds[b])
    min_pop <- if (is_null) bs else max(5L, floor(bs/5L))

    tr <- SimulateTrueCounts(
      ncells_total = bs, min_popsize = min_pop, i_minpop = 1L,
      ngenes = N_GENES, evf_type = evf_type,
      phyla = if (evf_type == "discrete") Phyla5() else NULL,
      randseed = bseed, Sigma = SYMSIM_SIGMA, n_de_evf = n_de_evf,
      bimod = 0, vary = "s"
    )
    all_counts[[b]] <- tr[["counts"]]
    all_pops[[b]]   <- as.character(tr[["cell_meta"]]$pop)
  }

  true_counts <- as(do.call(cbind, all_counts), "dgCMatrix")
  pop <- unlist(all_pops)
  rownames(true_counts) <- paste0("Gene", seq_len(nrow(true_counts)))

  # Full validation against already-verified production file (all 244, not sampled)
  prod_path <- file.path(PROD_DIR, sprintf("symsim_run_%05d.rds", run_id))
  matches <- NA
  if (file.exists(prod_path)) {
    prod <- readRDS(prod_path)
    matches <- identical(pop, prod$cell_meta$true_group)
    if (!isTRUE(matches)) {
      n_mismatch <- n_mismatch + 1L
      cat(sprintf("  *** MISMATCH vs production for run_id=%05d (fit_key=%s) ***\n", run_id, fk))
    }
  } else {
    cat(sprintf("  *** WARNING: production file not found: %s ***\n", prod_path))
  }

  group_levels <- sort(unique(pop))
  true_group_means <- vapply(group_levels, function(g)
    Matrix::rowMeans(true_counts[, pop==g, drop=FALSE]),
    FUN.VALUE = numeric(nrow(true_counts)))
  colnames(true_group_means) <- group_levels

  saveRDS(list(
    true_group_means = true_group_means,
    fit_key = fk,
    representative_run_id = run_id,
    pop_levels = group_levels,
    source = "symsim_true_counts_pre_observed",
    matches_production_pop = matches,
    method_note = paste("Group means of SimulateTrueCounts() output, before",
                         "True2ObservedCounts(). Validated by exact match against",
                         "already-verified production cell_meta$true_group.")
  ), out_path, compress=TRUE)

  manifest <- rbind(manifest, data.frame(fit_key=fk, representative_run_id=run_id,
                     n_populations=length(group_levels), n_cells_total=length(pop),
                     matches_production_pop=matches, output_path=out_path))

  cat(sprintf("[%3d/%d] DONE fit_key=%-30s run_id=%05d  pops=%d  match=%s\n",
              i, length(unique_fits), fk, run_id, length(group_levels), matches))
}

write.csv(manifest, "data/simulated/ground_truth/symsim_manifest.csv", row.names=FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
cat(sprintf("\n=== DONE: %d fit_keys in %.1f sec | mismatches: %d ===\n",
            length(unique_fits), elapsed, n_mismatch))
