# EDA 2.3 (closing check) — per-dataset PCA (log1p counts), colored by donor.
# Isolates within-technology batch structure without cross-technology depth
# confound (see combined raw/log1p results: dominated by a depth-driven
# "horseshoe" artifact, not directly informative about donor-level batch
# structure).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(irlba)
})

run_one <- function(name, path) {
  sce <- readRDS(path)
  cts <- counts(sce)
  log_cts <- cts
  log_cts@x <- log1p(log_cts@x)

  pca <- irlba(t(log_cts), nv = 10, center = Matrix::rowMeans(log_cts))
  var_explained <- pca$d^2 / sum(pca$d^2)
  cat(name, "- variance explained (first 3 PCs):", round(100*var_explained[1:3], 1), "\n")

  list(pca = pca, var_explained = var_explained, donor = as.character(sce$batch_id), name = name)
}

muraro <- run_one("muraro", "/mnt/extra/scrnaseq-data/real/muraro/muraro_harmonized.rds")
baron <- run_one("baron", "/mnt/extra/scrnaseq-data/real/baron/baron_harmonized.rds")
segerstolpe <- run_one("segerstolpe", "/mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_harmonized.rds")

png("results/eda_checkpoint2/pancreas_batch_pca_per_dataset.png", width = 1500, height = 500)
par(mfrow = c(1, 3))
for (res in list(muraro, baron, segerstolpe)) {
  donor_factor <- factor(res$donor)
  plot(res$pca$u[,1], res$pca$u[,2], col = rainbow(nlevels(donor_factor))[donor_factor], pch = 16, cex = 0.7,
       xlab = paste0("PC1 (", round(100*res$var_explained[1],1), "%)"),
       ylab = paste0("PC2 (", round(100*res$var_explained[2],1), "%)"),
       main = paste0(res$name, " (", nlevels(donor_factor), " donors)"))
}
dev.off()
cat("Saved: results/eda_checkpoint2/pancreas_batch_pca_per_dataset.png\n")
