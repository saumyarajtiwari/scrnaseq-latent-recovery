# Verification pass for Step 2.4 -- run as a SEPARATE process (memory safety).
suppressPackageStartupMessages(library(SingleCellExperiment))
check <- readRDS("/mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_sce.rds")
cat("Dimensions:", dim(check), "\n")
cat("Counts class:", class(assay(check,"counts"))[1], "\n")
cat("Any all-zero cells:", any(colSums(assay(check,"counts")) == 0), "\n")
cat("Unique cell_subclass:", length(unique(check$cell_subclass)), "\n")
cat("Unique cell_cluster:", length(unique(check$cell_cluster)), "\n")
