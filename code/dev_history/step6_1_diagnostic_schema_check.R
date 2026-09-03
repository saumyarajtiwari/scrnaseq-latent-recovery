# Step 6.1 diagnostic: confirm exact field names for batch/donor covariates
# and cell-type/label columns across all six real datasets, plus confirm
# column names in the Step 4 master results table and embedding manifest,
# before writing the Step 6.1 production script.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

real_files <- list(
  pbmc68k     = "data/real/pbmc68k/pbmc68k_harmonized.rds",
  muraro      = "data/real/muraro/muraro_harmonized.rds",
  baron       = "data/real/baron/baron_harmonized.rds",
  segerstolpe = "data/real/segerstolpe/segerstolpe_harmonized.rds",
  ts_lung     = "data/real/ts_lung/ts_lung_harmonized.rds",
  tasic2018   = "data/real/tasic2018/tasic2018_harmonized.rds"
)

for (nm in names(real_files)) {
  path <- real_files[[nm]]
  if (!file.exists(path)) {
    cat(sprintf("[%s] FILE NOT FOUND at %s\n", nm, path))
    next
  }
  sce <- readRDS(path)
  cat(sprintf("\n=== %s ===\n", nm))
  cat("colData columns: ", paste(colnames(colData(sce)), collapse = ", "), "\n")

  batch_candidates <- grep("batch|donor|individual|sample|patient|subject",
                            colnames(colData(sce)), ignore.case = TRUE, value = TRUE)
  cat("Candidate batch/donor columns: ",
      if (length(batch_candidates) > 0) paste(batch_candidates, collapse = ", ") else "NONE FOUND", "\n")
  for (cc in batch_candidates) {
    vals <- colData(sce)[[cc]]
    cat(sprintf("  %s -> %d unique levels: %s\n", cc, length(unique(vals)),
                paste(head(unique(vals), 10), collapse = ", ")))
  }

  ct_candidates <- grep("cell_type|true_group|label|annotation|cluster",
                         colnames(colData(sce)), ignore.case = TRUE, value = TRUE)
  cat("Candidate cell-type/label columns: ",
      if (length(ct_candidates) > 0) paste(ct_candidates, collapse = ", ") else "NONE FOUND", "\n")
  for (cc in ct_candidates) {
    cat(sprintf("  %s -> %d unique levels\n", cc, length(unique(colData(sce)[[cc]]))))
  }
  rm(sce); gc()
}

cat("\n=== step4_master_results_table.csv columns ===\n")
mrt <- read.csv("data/processed/step4_master_results_table.csv", nrows = 5)
print(colnames(mrt))

cat("\n=== embedding_manifest.csv columns ===\n")
manifest <- read.csv("data/processed/embedding_manifest.csv", nrows = 5)
print(colnames(manifest))

cat("\n=== real_data_inventory.csv columns ===\n")
inv <- read.csv("data/real_data_inventory.csv", nrows = 10)
print(colnames(inv))
print(inv)

cat("\n=== aricode package check ===\n")
if (requireNamespace("aricode", quietly = TRUE)) {
  cat("aricode already installed, version:", as.character(packageVersion("aricode")), "\n")
} else {
  cat("aricode NOT installed — will need install.packages('aricode')\n")
}
