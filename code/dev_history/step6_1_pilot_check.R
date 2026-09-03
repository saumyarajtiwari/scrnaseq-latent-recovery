suppressPackageStartupMessages({
  library(Matrix)
})

# --- Part A: inspect structure of one simulated + one real embedding file ---
manifest <- read.csv("data/processed/embedding_manifest.csv")

cat("=== Sample simulated embedding file (Raw PCA, Splatter run 371) ===\n")
sim_row <- manifest[manifest$source == "splatter" & manifest$run_id == 371 &
                     manifest$method == "raw_pca" & manifest$is_null_control == FALSE, ]
print(sim_row[, c("file_path","method","source","run_id")])
if (nrow(sim_row) > 0) {
  obj <- readRDS(sim_row$file_path[1])
  cat("Top-level names: ", paste(names(obj), collapse=", "), "\n")
  str(obj, max.level = 1)
}

cat("\n=== Sample real embedding file (Raw PCA, Baron) ===\n")
real_row <- manifest[manifest$source == "baron" & manifest$method == "raw_pca", ]
print(real_row[, c("file_path","method","source")])
if (nrow(real_row) > 0) {
  obj2 <- readRDS(real_row$file_path[1])
  cat("Top-level names: ", paste(names(obj2), collapse=", "), "\n")
  str(obj2, max.level = 1)
}

# --- Part B: k-means reproducibility pilot ---
if (!requireNamespace("aricode", quietly = TRUE)) {
  install.packages("aricode", repos = "https://cloud.r-project.org")
}
library(aricode)

master <- read.csv("data/processed/step4_master_results_table.csv")

pilot <- master[master$source == "splatter" & master$method == "raw_pca" &
                master$is_null_control == FALSE & !is.na(master$ari), ]
pilot <- head(pilot, 5)
print(pilot[, c("run_id","method","ari","n_groups")])

# We need true_group labels + embedding to test k-means reproduction.
# Placeholder — will fill in once Part A confirms where true_group and
# embedding actually live relative to each other.
cat("\nNOTE: k-means reproduction test deferred until Part A confirms file structure.\n")
