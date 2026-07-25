# EDA Checkpoint 2.3 — Batch Effect Visibility Check (pancreas datasets only)
# PCA on RAW (unnormalized) counts, as specified. NOTE: raw-count PCA across
# three different sequencing technologies (CEL-seq2, inDrop, Smart-seq2) with
# a ~64x range in mean depth (per EDA 2.1) is expected to be dominated by
# depth/technology differences rather than biological batch structure alone.
# This is documented explicitly, not treated as a defect.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

muraro <- readRDS("/mnt/extra/scrnaseq-data/real/muraro/muraro_harmonized.rds")
baron <- readRDS("/mnt/extra/scrnaseq-data/real/baron/baron_harmonized.rds")
segerstolpe <- readRDS("/mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_harmonized.rds")

rownames(muraro) <- sub("__chr.*$", "", rownames(muraro))
cat("Muraro rownames after suffix strip:", paste(head(rownames(muraro),3), collapse=", "), "\n")

common_genes <- Reduce(intersect, list(rownames(muraro), rownames(baron), rownames(segerstolpe)))
cat("Common genes across all 3:", length(common_genes), "\n")

m <- as(counts(muraro)[common_genes, ], "CsparseMatrix")
b <- as(counts(baron)[common_genes, ], "CsparseMatrix")
s <- as(counts(segerstolpe)[common_genes, ], "CsparseMatrix")

combined_counts <- cbind(m, b, s)
dataset_label <- c(rep("muraro", ncol(m)), rep("baron", ncol(b)), rep("segerstolpe", ncol(s)))
donor_label <- c(as.character(muraro$batch_id), as.character(baron$batch_id), as.character(segerstolpe$batch_id))
donor_label_full <- paste0(dataset_label, "_", donor_label)

cat("Combined dims:", dim(combined_counts), "\n")
cat("Cells per dataset:", table(dataset_label), "\n")

pca <- prcomp(t(as.matrix(combined_counts)), scale. = FALSE, center = TRUE, rank. = 10)
cat("PCA done. Variance explained (first 5 PCs):",
    round(100 * pca$sdev[1:5]^2 / sum(pca$sdev^2), 1), "\n")

saveRDS(list(pca = pca, dataset_label = dataset_label, donor_label_full = donor_label_full),
        "results/eda_checkpoint2/pancreas_batch_pca_result.rds")

png("results/eda_checkpoint2/pancreas_batch_pca.png", width = 1400, height = 700)
par(mfrow = c(1, 2))

dataset_colors <- c(muraro = "steelblue", baron = "darkorange", segerstolpe = "forestgreen")
plot(pca$x[,1], pca$x[,2], col = dataset_colors[dataset_label], pch = 16, cex = 0.6,
     xlab = paste0("PC1 (", round(100*pca$sdev[1]^2/sum(pca$sdev^2),1), "%)"),
     ylab = paste0("PC2 (", round(100*pca$sdev[2]^2/sum(pca$sdev^2),1), "%)"),
     main = "Colored by dataset/technology")
legend("topright", legend = names(dataset_colors), col = dataset_colors, pch = 16, bty = "n")

donor_factor <- factor(donor_label_full)
donor_colors <- rainbow(length(levels(donor_factor)))[donor_factor]
plot(pca$x[,1], pca$x[,2], col = donor_colors, pch = 16, cex = 0.6,
     xlab = paste0("PC1 (", round(100*pca$sdev[1]^2/sum(pca$sdev^2),1), "%)"),
     ylab = paste0("PC2 (", round(100*pca$sdev[2]^2/sum(pca$sdev^2),1), "%)"),
     main = paste0("Colored by donor (", length(levels(donor_factor)), " total batches)"))

dev.off()
cat("Saved: results/eda_checkpoint2/pancreas_batch_pca.png\n")
