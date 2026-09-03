# EDA 2.3 (supplementary) — log1p-transformed PCA, pancreas datasets combined.
# Uses irlba (truncated SVD on sparse matrix directly, no dense coercion)
# for speed/memory safety, unlike the raw-count prcomp() run.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

muraro <- readRDS("/mnt/extra/scrnaseq-data/real/muraro/muraro_harmonized.rds")
baron <- readRDS("/mnt/extra/scrnaseq-data/real/baron/baron_harmonized.rds")
segerstolpe <- readRDS("/mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_harmonized.rds")
rownames(muraro) <- sub("__chr.*$", "", rownames(muraro))

common_genes <- Reduce(intersect, list(rownames(muraro), rownames(baron), rownames(segerstolpe)))
m <- counts(muraro)[common_genes, ]; b <- counts(baron)[common_genes, ]; s <- counts(segerstolpe)[common_genes, ]
combined_counts <- cbind(m, b, s)
combined_log <- combined_counts
combined_log@x <- log1p(combined_log@x)

dataset_label <- c(rep("muraro", ncol(m)), rep("baron", ncol(b)), rep("segerstolpe", ncol(s)))
donor_label_full <- paste0(dataset_label, "_", c(as.character(muraro$batch_id), as.character(baron$batch_id), as.character(segerstolpe$batch_id)))

pca <- irlba(t(combined_log), nv = 10, center = Matrix::rowMeans(combined_log))
var_explained <- pca$d^2 / sum(pca$d^2)
cat("Variance explained (first 5 PCs):", round(100*var_explained[1:5], 1), "\n")

png("results/eda_checkpoint2/pancreas_batch_pca_log1p.png", width = 1400, height = 700)
par(mfrow = c(1, 2))
dataset_colors <- c(muraro="steelblue", baron="darkorange", segerstolpe="forestgreen")
plot(pca$u[,1], pca$u[,2], col = dataset_colors[dataset_label], pch = 16, cex = 0.6,
     xlab = paste0("PC1 (", round(100*var_explained[1],1), "%)"),
     ylab = paste0("PC2 (", round(100*var_explained[2],1), "%)"),
     main = "log1p counts: colored by dataset/technology")
legend("topright", legend = names(dataset_colors), col = dataset_colors, pch = 16, bty = "n")

donor_factor <- factor(donor_label_full)
plot(pca$u[,1], pca$u[,2], col = rainbow(nlevels(donor_factor))[donor_factor], pch = 16, cex = 0.6,
     xlab = paste0("PC1 (", round(100*var_explained[1],1), "%)"),
     ylab = paste0("PC2 (", round(100*var_explained[2],1), "%)"),
     main = paste0("log1p counts: colored by donor (", nlevels(donor_factor), " batches)"))
dev.off()
cat("Saved: results/eda_checkpoint2/pancreas_batch_pca_log1p.png\n")
