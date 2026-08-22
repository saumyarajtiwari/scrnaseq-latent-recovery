suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
  library(scater)
  library(igraph)
  library(irlba)
  library(Matrix)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(42)
LOG <- function(...) cat(format(Sys.time()), "-", ..., "\n")

LOG("Loading PBMC68k...")
sce <- readRDS("data/real/pbmc68k/pbmc68k_harmonized.rds")
LOG("Cells:", ncol(sce), "Genes:", nrow(sce))

t0 <- Sys.time()
LOG("Quick clustering for size-factor estimation (standard scran practice at this scale)...")
quick_clusters <- quickCluster(sce)
LOG("Quick clusters:", length(unique(quick_clusters)), "- elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")), 2), "min")

t1 <- Sys.time()
LOG("Computing pooling-based size factors...")
sce <- computeSumFactors(sce, clusters = quick_clusters)
LOG("Done - elapsed:", round(as.numeric(difftime(Sys.time(), t1, units="mins")), 2), "min")

t2 <- Sys.time()
LOG("Log-normalizing...")
sce <- logNormCounts(sce)
LOG("Done - elapsed:", round(as.numeric(difftime(Sys.time(), t2, units="mins")), 2), "min")

t3 <- Sys.time()
LOG("Modeling gene variance, selecting top 2000 HVGs...")
var_fit <- modelGeneVar(sce)
hvgs <- getTopHVGs(var_fit, n = 2000)
LOG("HVGs selected:", length(hvgs), "- elapsed:", round(as.numeric(difftime(Sys.time(), t3, units="mins")), 2), "min")

t4 <- Sys.time()
LOG("Computing PCA (irlba, 30 components, matching PBMC 3k's rank=30) on HVG subset...")
logcounts_hvg <- as.matrix(logcounts(sce)[hvgs, ])
pca_res <- irlba(t(logcounts_hvg), nv = 30, center = TRUE, scale = FALSE)
pca_embedding <- pca_res$u %*% diag(pca_res$d)
rownames(pca_embedding) <- colnames(sce)
reducedDim(sce, "PCA") <- pca_embedding
rm(logcounts_hvg); gc(verbose = FALSE)
LOG("PCA done - elapsed:", round(as.numeric(difftime(Sys.time(), t4, units="mins")), 2), "min")

t5 <- Sys.time()
LOG("Building SNN graph (k=20)...")
snn_graph <- buildSNNGraph(sce, k = 20, use.dimred = "PCA")
LOG("Graph built - elapsed:", round(as.numeric(difftime(Sys.time(), t5, units="mins")), 2), "min")

t6 <- Sys.time()
LOG("Detecting communities via Louvain (substituted for Walktrap: OOM-killed at 13.3GB RSS on this 65,690-node graph; Louvain (Blondel et al. 2008) is the field-standard scalable alternative, default in Seurat/Scanpy for datasets at this scale)...")
lv <- cluster_louvain(snn_graph)
clusters <- as.factor(membership(lv))
LOG("Clusters found:", length(unique(clusters)), "- elapsed:", round(as.numeric(difftime(Sys.time(), t6, units="mins")), 2), "min")
print(table(clusters))

# PBMC68k rownames are Ensembl gene IDs, not symbols (confirmed empty rowData,
# no symbol column available) - map canonical marker symbols to Ensembl IDs
# via org.Hs.eg.db before scoring (verified: all 9 markers map cleanly and
# all 9 Ensembl IDs are present in this dataset's rownames)
marker_symbols <- list(
  CD4_T    = c("IL7R", "CD3D"),
  CD8_T    = c("CD8A", "CD3D"),
  B_cell   = c("CD79A", "MS4A1"),
  NK       = c("GNLY", "NKG7"),
  Monocyte = c("CD14", "LYZ")
)
all_symbols <- unique(unlist(marker_symbols))
symbol_to_ensembl <- mapIds(org.Hs.eg.db, keys = all_symbols, column = "ENSEMBL", keytype = "SYMBOL")
marker_sets <- lapply(marker_symbols, function(syms) unname(symbol_to_ensembl[syms]))
LOG("Marker symbol-to-Ensembl mapping:")
print(symbol_to_ensembl)

LOG("Scoring clusters against canonical marker sets...")
logcounts_full <- logcounts(sce)
cluster_scores <- sapply(names(marker_sets), function(ct) {
  genes_present <- marker_sets[[ct]][marker_sets[[ct]] %in% rownames(sce)]
  if (length(genes_present) == 0) return(rep(NA_real_, length(levels(clusters))))
  expr <- as.matrix(logcounts_full[genes_present, , drop = FALSE])
  mean_expr_per_cell <- colMeans(expr)
  tapply(mean_expr_per_cell, clusters, mean)
})
rownames(cluster_scores) <- levels(clusters)
print(round(cluster_scores, 3))

# Preserve raw cluster membership for reproducibility / future refinement
sce$louvain_cluster <- clusters

# Hierarchical gating fix: a flat 5-way marker race conflates NK with cytotoxic
# CD8+ T cells, since both express GNLY/NKG7 (shared cytotoxic-effector
# markers) - discovered via inspection of this exact run's cluster_scores
# table (clusters with substantial CD3D-driven T-cell scores were still
# losing to NK due to elevated GNLY/NKG7). Standard fix (flow cytometry and
# scRNA-seq annotation literature): gate on CD3D first (T cell vs not),
# THEN subtype - NK's cytotoxic markers should never compete against a
# genuine CD3D+ T cell for cluster assignment.
cd3d_ensembl <- unname(mapIds(org.Hs.eg.db, keys="CD3D", column="ENSEMBL", keytype="SYMBOL"))
cd3d_expr <- as.numeric(logcounts_full[cd3d_ensembl, ])
cd3d_by_cluster <- tapply(cd3d_expr, clusters, mean)
LOG("Per-cluster mean CD3D expression (T-cell gate):")
print(round(cd3d_by_cluster, 3))

il7r_ensembl <- unname(mapIds(org.Hs.eg.db, keys="IL7R", column="ENSEMBL", keytype="SYMBOL"))
cd8a_ensembl <- unname(mapIds(org.Hs.eg.db, keys="CD8A", column="ENSEMBL", keytype="SYMBOL"))
il7r_by_cluster <- tapply(as.numeric(logcounts_full[il7r_ensembl, ]), clusters, mean)
cd8a_by_cluster <- tapply(as.numeric(logcounts_full[cd8a_ensembl, ]), clusters, mean)

is_t_cell <- cd3d_by_cluster > median(cd3d_by_cluster)
LOG("T-cell gate result per cluster (CD3D > median):")
print(is_t_cell)

cluster_to_celltype <- character(length(levels(clusters)))
names(cluster_to_celltype) <- levels(clusters)
for (cl in levels(clusters)) {
  if (is_t_cell[cl]) {
    cluster_to_celltype[cl] <- if (il7r_by_cluster[cl] >= cd8a_by_cluster[cl]) "CD4_T" else "CD8_T"
  } else {
    non_t_scores <- cluster_scores[cl, c("B_cell", "NK", "Monocyte")]
    cluster_to_celltype[cl] <- names(which.max(non_t_scores))
  }
}
LOG("Cluster-to-celltype assignment (hierarchical, post-CD3D-gate):")
print(cluster_to_celltype)

new_true_group <- unname(cluster_to_celltype[as.character(clusters)])
LOG("New label distribution:")
print(table(new_true_group, useNA = "always"))

sce$true_group <- new_true_group

OUT_PATH <- "data/real/pbmc68k/pbmc68k_harmonized_REANNOTATED.rds"
saveRDS(sce, OUT_PATH)
LOG("Saved to", OUT_PATH, "for verification (original file untouched)")
LOG("=== PBMC68k RE-ANNOTATION COMPLETE ===")
