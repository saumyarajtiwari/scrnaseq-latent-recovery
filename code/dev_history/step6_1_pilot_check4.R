suppressPackageStartupMessages(library(aricode))

manifest <- read.csv("data/processed/embedding_manifest.csv")

# Check: does batch_id ever contain NA/empty, and does AMI() have the same
# character-coercion issue as ARI() did?
row <- manifest[manifest$source == "baron" & manifest$method == "pca_raw", ]
obj <- readRDS(row$file_path[1])
tg <- obj$true_group
bid <- obj$batch_id
cat("batch_id class:", class(bid), " unique:", paste(unique(bid), collapse=","), "\n")

set.seed(42)
km <- kmeans(obj$embedding, centers = length(unique(tg)), nstart = 25)$cluster

# Test AMI with raw character batch_id (expect possible failure/coercion issue)
test1 <- tryCatch(aricode::AMI(km, bid), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("AMI with raw character batch_id:", test1, "\n")

# Test AMI with factor-coded batch_id (the fix pattern from ARI)
bid_int <- as.integer(factor(bid))
test2 <- tryCatch(aricode::AMI(km, bid_int), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("AMI with factor-coded batch_id:", test2, "\n")

# Check source_file path pattern for raw counts access
cat("\nsource_file for this embedding:", obj$source_file, "\n")
cat("Does source_file exist:", file.exists(obj$source_file), "\n")
if (file.exists(obj$source_file)) {
  src <- readRDS(obj$source_file)
  cat("class(src):", class(src), "\n")
  if (is(src, "SingleCellExperiment")) {
    cat("assayNames:", paste(assayNames(src), collapse=","), "\n")
    cts <- assay(src, "counts")
    cat("dim(counts):", paste(dim(cts), collapse=" x "), "\n")
    umi <- Matrix::colSums(cts)
    cat("UMI summary:", paste(round(summary(umi)), collapse=", "), "\n")
  }
}
