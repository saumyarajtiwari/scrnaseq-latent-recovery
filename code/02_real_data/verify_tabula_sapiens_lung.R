# Verification pass for Step 2.3 — run as a SEPARATE process from the
# download/save script (see memory note in download_tabula_sapiens_lung.R).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

check <- readRDS("~/Desktop/scrnaseq-latent-recovery/data/real/tabula_sapiens_lung/tabula_sapiens_lung_raw.rds")

cat("Dimensions:", dim(check), "\n")
cat("Counts class:", class(assay(check, "counts"))[1], "\n")
cat("Unique cell types:", length(unique(check$cell_type)), "\n")
cat("Unique donors:", length(unique(check$donor_id)), "\n")
cat("Unique assay values:", paste(unique(check$assay), collapse=", "), "\n")
cat("Any all-zero cells:", any(colSums(assay(check,"counts")) == 0), "\n")
