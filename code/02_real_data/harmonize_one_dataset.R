# Step 2.6 — Format Harmonization for a single real dataset
# Usage: Rscript harmonize_one_dataset.R <in_path> <out_path> <dataset_name> \
#          <true_group_col_or_NA> <batch_col_or_NA> <extra_cols_comma_or_NA>
#
# Target schema (matched to Step 1 simulated SCEs where applicable):
#   assay: counts only (raw, unnormalized -- normalization happens in Step 3)
#   colData: cell_id, true_group, batch_id, achieved_sparsity, dataset_source,
#            plus dataset-specific extra columns preserved under their
#            original names (not force-mapped, since they have no simulated
#            equivalent: e.g. disease, cell_class, dissected_region)
#   metadata(): dataset_name, ground_truth_source, simulator (NA for real
#            data), harmonization_notes
#
# true_group is NA for pbmc68k (no annotation exists yet -- Step 2.1 decision,
# re-annotation deferred). This is NOT a bug; it's an explicit, documented gap.

args <- commandArgs(trailingOnly = TRUE)
in_path <- args[1]; out_path <- args[2]; dataset_name <- args[3]
true_group_col <- args[4]; batch_col <- args[5]; extra_cols <- args[6]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

sce <- readRDS(in_path)
cat("Loaded", dataset_name, "| dims:", dim(sce), "\n")

cts <- assay(sce, "counts")
if (!is(cts, "dgCMatrix")) cts <- as(cts, "CsparseMatrix")

achieved_sparsity <- 1 - (Matrix::nnzero(cts) / (as.numeric(nrow(cts)) * as.numeric(ncol(cts))))

true_group <- if (true_group_col == "NA") {
  factor(rep(NA, ncol(sce)))
} else {
  as.character(sce[[true_group_col]])
}

batch_id <- if (batch_col == "NA") {
  factor(rep(NA, ncol(sce)))
} else {
  as.character(sce[[batch_col]])
}

cell_id_vec <- colnames(sce)
if (length(cell_id_vec) == 0) {
  cat("colnames(sce) empty - falling back to Barcode column for cell_id
")
  cell_id_vec <- sce$Barcode
}
stopifnot(length(cell_id_vec) == ncol(sce))

new_coldata <- DataFrame(
  cell_id = cell_id_vec,
  true_group = true_group,
  batch_id = batch_id,
  achieved_sparsity = achieved_sparsity,
  dataset_source = "real"
)

if (extra_cols != "NA") {
  for (col in strsplit(extra_cols, ",")[[1]]) {
    new_coldata[[col]] <- sce[[col]]
  }
}

ground_truth_note <- if (true_group_col == "NA") {
  "No published cell-type annotation available; re-annotation deferred (Step 2.1 decision)."
} else {
  paste0("Published author cell-type annotation, original column: '", true_group_col, "'")
}

sce_harmonized <- SingleCellExperiment(
  assays = list(counts = cts),
  colData = new_coldata
)
rownames(sce_harmonized) <- rownames(sce)

metadata(sce_harmonized) <- list(
  dataset_name = dataset_name,
  ground_truth_source = ground_truth_note,
  simulator = NA,
  harmonization_notes = paste0(
    "Harmonized Step 2.6. Original batch column: '", batch_col,
    "'. Extra preserved columns: '", extra_cols, "'."
  )
)

cat("Harmonized dims:", dim(sce_harmonized), "\n")
cat("Assay names:", assayNames(sce_harmonized), "\n")
cat("colData columns:", paste(colnames(colData(sce_harmonized)), collapse=", "), "\n")

saveRDS(sce_harmonized, out_path, compress = TRUE)
cat("Saved:", out_path, "\n")
