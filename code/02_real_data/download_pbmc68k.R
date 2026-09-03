# Step 2.1 — PBMC 68k Download
# Downloads the Fresh 68k PBMCs (Donor A) dataset via TENxPBMCData,
# materializes the DelayedMatrix-backed counts to an in-memory dgCMatrix
# (required before saving — see Step 1.4 precedent for the same issue),
# and saves to the project's real-data directory.
#
# Note: this dataset ships with NO manually-curated or bundled cell-type
# annotation (colData contains only technical/sample fields). The official
# "annotation" available elsewhere for this dataset is a classifier-inferred
# soft label against 10 FACS-purified references, with a companion k-means
# cluster file covering only 40k/68k cells with inconsistent numbering.
# Decision (this project): re-annotate independently later using the same
# marker-based approach applied to PBMC 3k in Step 1.4, for methodological
# consistency and defensibility. Not done in this script.

library(TENxPBMCData)
library(SingleCellExperiment)
library(Matrix)

pbmc68k <- TENxPBMCData(dataset = "pbmc68k")
counts(pbmc68k) <- as(counts(pbmc68k), "dgCMatrix")

out_dir <- "~/Desktop/scrnaseq-latent-recovery/data/real/pbmc68k"
saveRDS(pbmc68k, file.path(out_dir, "pbmc68k_raw.rds"), compress = TRUE)

cat("Saved PBMC 68k raw counts.\n")
cat("Dimensions:", dim(pbmc68k), "\n")
cat("Sparsity:", 1 - (nnzero(counts(pbmc68k)) / (as.numeric(nrow(pbmc68k)) * as.numeric(ncol(pbmc68k)))), "\n")
