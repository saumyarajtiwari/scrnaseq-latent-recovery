# Step 2.8 — Final Dataset Inventory
# Builds a machine-readable master inventory table, one row per real dataset,
# extracted directly from the harmonized SCE objects (not transcribed from
# the Step 2.7 catalog, to avoid two sources of truth drifting apart).
# This table is intended to be read programmatically by downstream pipeline
# scripts (Step 3 onward), not just as human documentation.

suppressPackageStartupMessages(library(SingleCellExperiment))

paths <- list(
  pbmc68k = "/mnt/extra/scrnaseq-data/real/pbmc68k/pbmc68k_harmonized.rds",
  muraro = "/mnt/extra/scrnaseq-data/real/muraro/muraro_harmonized.rds",
  baron = "/mnt/extra/scrnaseq-data/real/baron/baron_harmonized.rds",
  segerstolpe = "/mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_harmonized.rds",
  ts_lung = "/mnt/extra/scrnaseq-data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds",
  tasic2018 = "/mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_harmonized.rds"
)

rows <- list()
for (name in names(paths)) {
  path <- paths[[name]]
  sce <- readRDS(path)

  n_cell_types <- length(unique(sce$true_group[!is.na(sce$true_group)]))
  n_unlabeled <- sum(is.na(sce$true_group))
  fully_labeled <- (n_unlabeled == 0)

  rows[[name]] <- data.frame(
    dataset_name = name,
    n_cells = ncol(sce),
    n_genes = nrow(sce),
    n_cell_types = n_cell_types,
    n_cells_unlabeled = n_unlabeled,
    fully_labeled = fully_labeled,
    n_batches = length(unique(sce$batch_id)),
    approx_sparsity = round(mean(sce$achieved_sparsity), 4),
    file_path = path,
    stringsAsFactors = FALSE
  )
}

inventory <- do.call(rbind, rows)
rownames(inventory) <- NULL

cat("=== Final Inventory ===\n")
print(inventory)

out_path <- "~/Desktop/scrnaseq-latent-recovery/data/real_data_inventory.csv"
write.csv(inventory, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
