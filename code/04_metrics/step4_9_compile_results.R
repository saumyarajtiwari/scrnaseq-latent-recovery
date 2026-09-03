suppressPackageStartupMessages({})

primary <- read.csv("data/processed/step4_1to3_subspace_metrics.csv", stringsAsFactors = FALSE)
secondary <- read.csv("data/processed/step4_4to6_secondary_metrics.csv", stringsAsFactors = FALSE)
manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)

names(primary)[names(primary) == "status"] <- "status_primary"
names(primary)[names(primary) == "error_msg"] <- "error_msg_primary"
names(secondary)[names(secondary) == "status"] <- "status_secondary"
names(secondary)[names(secondary) == "error_msg"] <- "error_msg_secondary"

# secondary is the superset (covers both simulated + real); left-join primary
# onto it, so real rows naturally get NA for primary (geometric) metrics,
# consistent with Step 4.8's spec (no ground-truth subspace for real data)
combined <- merge(secondary, primary, by = c("source", "run_id", "method"), all.x = TRUE)

# Bring in parameter-grid / real-inventory metadata from the manifest,
# excluding file-integrity-scan-specific columns not relevant to the
# metrics table itself
manifest_meta_cols <- c("data_type", "source", "run_id", "method",
                         "sparsity", "depth", "dropout", "separability", "n_cells", "n_groups",
                         "batch", "gene_strategy", "clipping", "is_null_control",
                         "n_genes", "n_cell_types", "n_cells_unlabeled", "fully_labeled",
                         "n_batches", "approx_sparsity")
manifest_meta <- unique(manifest[, manifest_meta_cols])

final <- merge(combined, manifest_meta, by = c("data_type", "source", "run_id", "method"), all.x = TRUE)

OUT_PATH <- "data/processed/step4_master_results_table.csv"
write.csv(final, OUT_PATH, row.names = FALSE)

cat("=== STEP 4.9 COMPILATION COMPLETE ===\n")
cat("Total rows:", nrow(final), "\n")
cat("Total columns:", ncol(final), "\n")
cat("\nBy data_type:\n"); print(table(final$data_type))
cat("\nSimulated rows with primary metrics present:", sum(final$data_type=="simulated" & !is.na(final$grassmann_distance)), "\n")
cat("Real rows with primary metrics (should be 0, no ground truth):", sum(final$data_type=="real" & !is.na(final$grassmann_distance)), "\n")
cat("\nColumn names:\n"); print(names(final))
