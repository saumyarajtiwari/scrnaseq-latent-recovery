# =============================================================================
# calibrate_scdesign3.R
# Phase A calibration for scDesign3 simulation runs.
#
# Tasks:
#   1. Load and QC PBMC 3k reference
#   2. Normalize and select top 2000 HVGs
#   3. Cluster and annotate 5 cell types with canonical markers
#   4. Verify separability subsets have sufficient cells
#   5. Verify depth multipliers produce expected UMI/cell levels
#   6. Save annotated reference for simulate_scdesign3.R
#
# Output:
#   data/simulated/pbmc3k_annotated.rds
#   data/simulated/scdesign3_calib_depth.csv
# =============================================================================

suppressPackageStartupMessages({
  library(TENxPBMCData)
  library(scran)
  library(scuttle)
  library(bluster)
  library(BiocSingular)
  library(Matrix)
  library(SingleCellExperiment)
})

set.seed(42)

cat("=== scDesign3 Calibration ===\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))

# =============================================================================
# 1. LOAD REFERENCE
# =============================================================================

cat("1. Loading PBMC 3k reference...\n")
pbmc <- TENxPBMCData('pbmc3k')
cat(sprintf("   Raw: %d genes x %d cells\n", nrow(pbmc), ncol(pbmc)))

# =============================================================================
# 2. QC FILTERING
# =============================================================================

cat("\n2. QC filtering...\n")

pbmc <- addPerCellQCMetrics(pbmc)

# Standard PBMC 3k thresholds
keep_cells <- pbmc$detected >= 200 &
              pbmc$detected <= 2500 &
              pbmc$sum      >= 500
pbmc <- pbmc[, keep_cells]
cat(sprintf("   After QC: %d cells retained\n", ncol(pbmc)))

# =============================================================================
# 3. NORMALIZATION AND HVG SELECTION
# =============================================================================

cat("\n3. Normalization and HVG selection...\n")

clusters_quick <- quickCluster(pbmc, BSPARAM=IrlbaParam())
pbmc <- computeSumFactors(pbmc, clusters=clusters_quick)
pbmc <- logNormCounts(pbmc)

hvg_stats <- modelGeneVar(pbmc)
hvgs      <- getTopHVGs(hvg_stats, n=2000)
cat(sprintf("   Selected %d HVGs\n", length(hvgs)))

pbmc_hvg <- pbmc[hvgs, ]

# =============================================================================
# 4. PCA AND CLUSTERING
# =============================================================================

cat("\n4. PCA and clustering...\n")

mat      <- t(as.matrix(logcounts(pbmc_hvg)))
pca_res  <- prcomp(mat, rank.=30, scale.=FALSE)
reducedDim(pbmc_hvg, 'PCA') <- pca_res$x

g        <- buildSNNGraph(pbmc_hvg, use.dimred='PCA', k=20)
clusters <- igraph::cluster_walktrap(g)$membership
pbmc_hvg$cluster <- factor(clusters)

cat(sprintf("   Found %d clusters\n", length(unique(clusters))))
cat("   Cluster sizes:\n")
print(table(pbmc_hvg$cluster))

# =============================================================================
# 5. CELL TYPE ANNOTATION
# =============================================================================

cat("\n5. Annotating cell types with canonical markers...\n")

# Set gene symbol rownames
gene_symbols <- rowData(pbmc_hvg)$Symbol_TENx
if (is.null(gene_symbols) || all(is.na(gene_symbols))) {
  gene_symbols <- rowData(pbmc_hvg)$Symbol
}
rownames(pbmc_hvg) <- make.unique(as.character(gene_symbols))

markers <- list(
  "CD4_T"    = c("IL7R",  "CD3D"),
  "CD8_T"    = c("CD8A",  "CD3D"),
  "B_cell"   = c("CD79A", "MS4A1"),
  "NK"       = c("GNLY",  "NKG7"),
  "Monocyte" = c("CD14",  "LYZ")
)

cat("   Marker gene presence in HVGs:\n")
for (ct in names(markers)) {
  present <- markers[[ct]][markers[[ct]] %in% rownames(pbmc_hvg)]
  cat(sprintf("   %-12s: %s\n", ct,
      if (length(present)) paste(present, collapse=", ") else "NONE FOUND"))
}

# Score each cluster for each cell type
log_counts  <- logcounts(pbmc_hvg)
cluster_ids <- levels(pbmc_hvg$cluster)

score_mat <- matrix(
  0,
  nrow = length(cluster_ids),
  ncol = length(markers),
  dimnames = list(cluster_ids, names(markers))
)

for (ct in names(markers)) {
  genes_present <- markers[[ct]][markers[[ct]] %in% rownames(pbmc_hvg)]
  if (length(genes_present) > 0) {
    for (cl in cluster_ids) {
      cells_cl <- pbmc_hvg$cluster == cl
      score_mat[cl, ct] <- mean(
        colMeans(log_counts[genes_present, cells_cl, drop=FALSE])
      )
    }
  }
}

cat("\n   Marker scores per cluster (higher = stronger signal):\n")
print(round(score_mat, 3))

# Greedy unique assignment: each cell type assigned to best-matching cluster first
# Remaining clusters (when n_clusters > n_cell_types) get next-best available
score_df   <- as.data.frame(score_mat)
remaining  <- rownames(score_df)
ct_pool    <- colnames(score_df)
assignment <- character(nrow(score_df))
names(assignment) <- rownames(score_df)

for (round in seq_along(ct_pool)) {
  if (length(ct_pool) == 0 || length(remaining) == 0) break
  sub      <- score_df[remaining, ct_pool, drop=FALSE]
  idx      <- which(sub == max(sub), arr.ind=TRUE)[1, ]
  best_cl  <- rownames(sub)[idx[1]]
  best_ct  <- colnames(sub)[idx[2]]
  assignment[best_cl] <- best_ct
  remaining <- remaining[remaining != best_cl]
  ct_pool   <- ct_pool[ct_pool   != best_ct]
}
for (cl in remaining) {
  assignment[cl] <- names(which.max(score_mat[cl, ]))
}

cluster_to_celltype <- assignment
pbmc_hvg$cell_type  <- cluster_to_celltype[as.character(pbmc_hvg$cluster)]

cat("\n   Cluster -> cell type assignments:\n")
for (cl in cluster_ids) {
  cat(sprintf("   Cluster %-3s -> %s\n", cl, cluster_to_celltype[cl]))
}

cat("\n   Final cell type counts:\n")
print(table(pbmc_hvg$cell_type))

# =============================================================================
# 6. VERIFY SEPARABILITY SUBSETS
# =============================================================================

cat("\n6. Verifying separability subsets...\n")

separability_subsets <- list(
  "null"   = c("CD4_T"),
  "low"    = c("CD4_T", "CD8_T"),
  "medium" = c("CD4_T", "CD8_T", "B_cell", "NK"),
  "high"   = c("CD4_T", "CD8_T", "B_cell", "NK", "Monocyte")
)

all_types_present <- TRUE
for (level in names(separability_subsets)) {
  types     <- separability_subsets[[level]]
  cells_in  <- pbmc_hvg$cell_type %in% types
  n_cells   <- sum(cells_in)
  n_types   <- length(unique(pbmc_hvg$cell_type[cells_in]))
  min_count <- min(table(pbmc_hvg$cell_type[cells_in]))
  flag      <- if (min_count < 50) " << WARNING: low cell count" else ""
  cat(sprintf("   %-8s %d cells, %d type(s), min per type = %d%s\n",
              level, n_cells, n_types, min_count, flag))
  if (n_types != length(types)) all_types_present <- FALSE
}

if (!all_types_present) {
  cat("\n   WARNING: Some expected cell types were not found.\n")
  cat("   Review cluster assignments above and adjust annotation if needed.\n")
} else {
  cat("\n   All separability subsets OK.\n")
}

# =============================================================================
# 7. VERIFY DEPTH MULTIPLIERS
# =============================================================================

cat("\n7. Verifying depth multipliers...\n")

baseline_depth <- mean(colSums(counts(pbmc_hvg)))
cat(sprintf("   Baseline mean depth (2000 HVGs): %.1f UMI/cell\n\n", baseline_depth))

targets    <- c(500, 2000, 10000)
multipliers <- round(targets / baseline_depth, 4)
depth_df <- data.frame(
  depth_label    = c("500", "2000", "10000"),
  multiplier     = multipliers,
  expected_depth = round(baseline_depth * multipliers, 1)
)

for (i in seq_len(nrow(depth_df))) {
  cat(sprintf("   multiplier %.2f -> expected %.1f UMI/cell  (target: %s)\n",
              depth_df$multiplier[i],
              depth_df$expected_depth[i],
              depth_df$depth_label[i]))
}

# =============================================================================
# 8. SAVE OUTPUTS
# =============================================================================

cat("\n8. Saving outputs...\n")

saveRDS(pbmc_hvg, "data/simulated/pbmc3k_annotated.rds")
cat("   Saved: data/simulated/pbmc3k_annotated.rds\n")

write.csv(depth_df, "data/simulated/scdesign3_calib_depth.csv", row.names=FALSE)
cat("   Saved: data/simulated/scdesign3_calib_depth.csv\n")

cat(sprintf("\n=== Calibration complete: %s ===\n", format(Sys.time())))
cat("Review cell type assignments and depth multipliers above\n")
cat("before proceeding to simulate_scdesign3.R\n")
