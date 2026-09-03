suppressPackageStartupMessages(library(SingleCellExperiment))

inv <- read.csv("data/real_data_inventory.csv")
print(inv[, c("dataset_name", "file_path", "n_batches")])

ts_path <- inv$file_path[inv$dataset_name == "ts_lung"]
cat("\nChecking ts_lung at:", ts_path, "\n")
if (file.exists(ts_path)) {
  sce <- readRDS(ts_path)
  cat("colData columns: ", paste(colnames(colData(sce)), collapse = ", "), "\n")
  bcol <- grep("batch|donor|individual|sample|patient|subject",
               colnames(colData(sce)), ignore.case = TRUE, value = TRUE)
  for (cc in bcol) {
    vals <- colData(sce)[[cc]]
    cat(sprintf("  %s -> %d unique levels: %s\n", cc, length(unique(vals)),
                paste(head(unique(vals), 10), collapse = ", ")))
  }
  rm(sce); gc()
} else {
  cat("STILL NOT FOUND at that path either.\n")
}

cat("\n=== Archive SSD mount check ===\n")
system("mount | grep archive")
system("df -h /mnt/archive 2>&1")

cat("\n=== Confirm libnorm/log/shiftedlog symlinks resolve ===\n")
for (m in c("pca_libnorm", "pca_log", "pca_shiftedlog")) {
  p <- file.path("data/processed", m)
  cat(m, "->", Sys.readlink(p), " exists:", dir.exists(p), "\n")
}
