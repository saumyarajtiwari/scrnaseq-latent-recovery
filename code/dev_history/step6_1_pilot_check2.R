manifest <- read.csv("data/processed/embedding_manifest.csv")
master   <- read.csv("data/processed/step4_master_results_table.csv")

cat("=== embedding_manifest.csv: unique values ===\n")
cat("method: ", paste(sort(unique(manifest$method)), collapse=", "), "\n")
cat("source: ", paste(sort(unique(manifest$source)), collapse=", "), "\n")
cat("data_type: ", paste(sort(unique(manifest$data_type)), collapse=", "), "\n")
cat("is_null_control: ", paste(sort(unique(manifest$is_null_control)), collapse=", "), "\n")
cat("status: ", paste(sort(unique(manifest$status)), collapse=", "), "\n")

cat("\n=== step4_master_results_table.csv: unique values ===\n")
cat("method: ", paste(sort(unique(master$method)), collapse=", "), "\n")
cat("source: ", paste(sort(unique(master$source)), collapse=", "), "\n")
cat("data_type: ", paste(sort(unique(master$data_type)), collapse=", "), "\n")
cat("is_null_control: ", paste(sort(unique(master$is_null_control)), collapse=", "), "\n")

cat("\n=== First 3 rows of each, all columns ===\n")
print(head(manifest, 3))
print(head(master, 3))
