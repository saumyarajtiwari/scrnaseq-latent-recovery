# Step 2.3 — Tabula Sapiens Lung Subset Download
# Downloads the Tabula Sapiens Lung dataset specifically (not all lung data
# in the Census, which spans many unrelated studies) via CZ CELLxGENE Census,
# using the Homo sapiens Census pinned to version 2025-11-08 for reproducibility.
#
# dataset_id "0d2ee4ac-05ee-40b2-afb6-ebb584caa867" identified via metadata
# query as the unique Tabula Sapiens Lung entry (65,847 cells, confirmed
# matching the Census-wide dataset table before any cell data was pulled).
#
# This subset was selected (per project scope) for its broad, sparse spread
# of 34 distinct cell types across only 4 donors, rather than a small number
# of dense, well-separated clusters. Note it also contains an internal
# technology batch axis (Smart-seq2 and 10x 3' v3 within the same dataset_id),
# preserved here via the `assay` column, in addition to donor_id.
#
# MEMORY NOTE: this dataset (~645MB saved, ~10GB+ RAM during conversion) has
# previously pushed this machine to ~94% RAM + swap usage. Do NOT hold the
# pulled object and a re-read verification copy in memory simultaneously in
# one process. This script deliberately runs pull+save and verification as
# two separate Rscript invocations so each fully exits and releases memory
# before the next begins.

suppressPackageStartupMessages({
  library(cellxgene.census)
  library(SingleCellExperiment)
  library(Matrix)
})

census <- open_soma(census_version = "2025-11-08")

sce <- get_single_cell_experiment(
  census = census,
  organism = "Homo sapiens",
  obs_value_filter = "dataset_id == '0d2ee4ac-05ee-40b2-afb6-ebb584caa867'",
  obs_column_names = c("cell_type", "tissue", "tissue_general", "assay",
                       "donor_id", "disease", "sex", "suspension_type")
)
census$close()

counts(sce) <- as(counts(sce), "CsparseMatrix")

out_path <- "~/Desktop/scrnaseq-latent-recovery/data/real/tabula_sapiens_lung/tabula_sapiens_lung_raw.rds"
saveRDS(sce, out_path, compress = TRUE)
cat("Saved:", out_path, "\n")
cat("Dimensions:", dim(sce), "\n")
