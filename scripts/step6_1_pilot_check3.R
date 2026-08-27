manifest <- read.csv("data/processed/embedding_manifest.csv")

cat("=== Sample simulated embedding file (Raw PCA, Splatter, run 371) ===\n")
sim_row <- manifest[manifest$source == "splatter" & manifest$method == "pca_raw" &
                     manifest$run_id == 371 & manifest$is_null_control == FALSE, ]
print(sim_row[, c("file_path","method","source","run_id")])
if (nrow(sim_row) > 0) {
  obj <- readRDS(sim_row$file_path[1])
  cat("Top-level names: ", paste(names(obj), collapse=", "), "\n")
  str(obj, max.level = 2, list.len = 20)
}

cat("\n=== Sample real embedding file (Raw PCA, Baron) ===\n")
real_row <- manifest[manifest$source == "baron" & manifest$method == "pca_raw", ]
print(real_row[, c("file_path","method","source")])
if (nrow(real_row) > 0) {
  obj2 <- readRDS(real_row$file_path[1])
  cat("Top-level names: ", paste(names(obj2), collapse=", "), "\n")
  str(obj2, max.level = 2, list.len = 20)
}
