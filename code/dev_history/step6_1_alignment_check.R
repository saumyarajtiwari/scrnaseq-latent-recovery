suppressPackageStartupMessages(library(SingleCellExperiment))

manifest <- read.csv("data/processed/embedding_manifest.csv")
row <- manifest[manifest$source == "baron" & manifest$method == "pca_raw", ]
obj <- readRDS(row$file_path[1])
src <- readRDS(obj$source_file)

cat("n cells in embedding object's true_group:", length(obj$true_group), "\n")
cat("n cells in source_file (ncol):", ncol(src), "\n")

# Compare batch_id order directly, if source has it
src_batch <- colData(src)$batch_id
cat("source colData has batch_id:", !is.null(src_batch), "\n")
if (!is.null(src_batch)) {
  cat("identical batch_id order:", identical(as.character(src_batch), as.character(obj$batch_id)), "\n")
  cat("n mismatches:", sum(as.character(src_batch) != as.character(obj$batch_id)), "\n")
}

# Compare true_group order too
src_tg <- colData(src)$true_group
if (!is.null(src_tg)) {
  cat("identical true_group order:", identical(as.character(src_tg), as.character(obj$true_group)), "\n")
}

# Also check cell_id / colnames alignment if present
cat("colnames(src) head:", paste(head(colnames(src), 3), collapse=","), "\n")
cat("colData(src)$cell_id head:", paste(head(colData(src)$cell_id, 3), collapse=","), "\n")
