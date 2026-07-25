# Step 2.5 — QC filter for a single real dataset (invoked separately per dataset for memory safety)
# Usage: Rscript qc_filter_one_dataset.R <in_path> <out_path> <mito_mode> <block_col_or_NA>
# mito_mode: "ensembl_pbmc68k" | "ensembl_tslung" | "none"

args <- commandArgs(trailingOnly = TRUE)
in_path <- args[1]; out_path <- args[2]; mito_mode <- args[3]; block_col <- args[4]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
})

sce <- readRDS(in_path)
cat("Loaded:", in_path, "| dims:", dim(sce), "\n")

human_mt_ensembl <- c("ENSG00000198888","ENSG00000198763","ENSG00000198804",
                      "ENSG00000198712","ENSG00000228253","ENSG00000198899",
                      "ENSG00000198938","ENSG00000198840","ENSG00000212907",
                      "ENSG00000198886","ENSG00000198786","ENSG00000198695","ENSG00000198727")

is_mito <- switch(mito_mode,
  "ensembl_pbmc68k" = rowData(sce)$ENSEMBL_ID %in% human_mt_ensembl,
  "ensembl_tslung" = rownames(sce) %in% human_mt_ensembl,
  "none" = NULL
)

subsets <- if (is.null(is_mito)) list() else list(Mito = is_mito)
qc <- perCellQCMetrics(sce, subsets = subsets)

batch_arg <- if (block_col == "NA") NULL else sce[[block_col]]
if (!is.null(batch_arg)) cat("Blocking QC by:", block_col, "(", length(unique(batch_arg)), "groups)\n")

low_sum <- isOutlier(qc$sum, log = TRUE, type = "lower", nmads = 3, batch = batch_arg)
low_detected <- isOutlier(qc$detected, log = TRUE, type = "lower", nmads = 3, batch = batch_arg)
cell_flagged <- low_sum | low_detected
if (!is.null(is_mito)) {
  high_mito <- isOutlier(qc$subsets_Mito_percent, type = "higher", nmads = 3, batch = batch_arg)
  cell_flagged <- cell_flagged | high_mito
}
cat("Cells flagged:", sum(cell_flagged), "of", ncol(sce), "(", round(100*sum(cell_flagged)/ncol(sce),1), "%)\n")

gene_detected <- rowSums(counts(sce) > 0)
gene_flagged <- isOutlier(gene_detected, log = TRUE, type = "lower", nmads = 3)
cat("Genes flagged:", sum(gene_flagged), "of", nrow(sce), "\n")

sce_filtered <- sce[!gene_flagged, !cell_flagged]
cat("Filtered dims:", dim(sce_filtered), "\n")

saveRDS(sce_filtered, out_path, compress = TRUE)
cat("Saved:", out_path, "\n")
