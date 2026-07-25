# Step 2.4 — Brain Atlas Subset Download (Tasic et al. 2018, GSE115746)
# Mouse visual cortex (VISp) + anterior lateral motor cortex (ALM), Smart-seq2.
#
# UNLIKE Steps 2.1-2.3, this dataset has NO Bioconductor/ExperimentHub wrapper.
# NOTE: scRNAseq::TasicBrainData() is a DIFFERENT, EARLIER dataset (Tasic et al.
# 2016, GSE71585, ~1800 cells) -- do not confuse the two.
#
# Source files (direct GEO FTP download, resumable):
#   GSE115746_cells_exon_counts.csv.gz         (45768 genes x 23178 cells)
#   GSE115746_complete_metadata_28706-cells.csv.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/
#
# DOWNLOAD LESSON: GEO's HTTPS download endpoint does not support byte-range
# resume (curl error 33). The FTP mirror does support resume, but an initial
# resumed-FTP-transfer attempt produced a file that matched the expected byte
# count exactly yet FAILED `gzip -t` integrity check (likely an ASCII/binary
# transfer-mode mismatch across the resume boundary). Fix: single-shot,
# non-resumed download with --ftp-method nocwd. ALWAYS run `gzip -t` on any
# FTP-downloaded .gz before trusting it, regardless of whether the byte count
# matches the expected size.
#
# MEMORY LESSON: a naive single-pass `fread` + `as.matrix` + sparsify of the
# full 45768 x 23178 exon-counts file (~1.06 billion cells, ~4.2GB as raw
# integers) was OOM-killed on this 15GB-RAM machine (Vcells peaked ~12.65GB
# before the kill). Fixed by reading and sparsifying in chunks of 5000 genes
# at a time (each chunk ~463MB dense), rbind-ing the resulting sparse pieces.
# Peak memory with chunking: ~7.8GB. Do not revert to a single-pass read.
#
# CELL-COUNT / CLASS-COUNT DISCREPANCY (investigated, not fixed further):
# Original scoping in this project referenced "23,822 cells, 133 types,
# 11 top-level classes." None of these three numbers exactly match what
# exists in the actual downloaded, published data:
#   - Real Mus musculus cells excluding controls/spike-ins/Low-Quality/
#     No-Class (per `cell_class` and `organism` columns in the metadata):
#     24,061 cells.
#   - Of those 24,061, only 22,113 have a matching column in the exon-counts
#     matrix. The other 1,948 are spread proportionally across nearly every
#     cell_subclass and both dissected_region values, and ALL pass QC
#     (sequencing_qc_pass_fail == "Pass") -- ruling out an ID-matching bug
#     (0 duplicate IDs on either side) or a systematic exclusion. Most
#     plausible explanation: these cells were profiled and passed QC, but
#     their exon-level counts were not included in this specific GEO
#     supplementary file (possibly intron-only, or a later-added cohort).
#     Not independently confirmed beyond this pattern-based inference.
#   - No metadata column contains exactly 11 categories. `cell_class` has
#     4 (Endothelial, GABAergic, Glutamatergic, Non-Neuronal). `cell_subclass`
#     has 24. The source of the "11 classes" figure could not be traced to
#     anything in the actual downloaded data; decision made (by delegated
#     project judgment) to use `cell_subclass` (24 categories) as the
#     top-level grouping for primary analysis, since it is an official,
#     citable Allen Institute metadata column, rather than construct an
#     unconfirmed custom 11-way merge of subclasses.
#   - `cell_cluster` has 134 unique values (vs. spec's 133) -- not
#     individually investigated; consistent with an unclassified catch-all
#     category appearing elsewhere in this same metadata file.
#
# FINAL DATASET: 45,768 genes x 22,113 cells. Top-level class = cell_subclass
# (24 categories). Fine-grained type = cell_cluster (134 categories).

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
})

raw_dir <- "/mnt/extra/scrnaseq-data/real/tasic2018/raw_download"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

exon_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746%5Fcells%5Fexon%5Fcounts%2Ecsv%2Egz"
meta_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115746&format=file&file=GSE115746%5Fcomplete%5Fmetadata%5F28706%2Dcells%2Ecsv%2Egz"

exon_path <- file.path(raw_dir, "GSE115746_cells_exon_counts.csv.gz")
meta_path <- file.path(raw_dir, "GSE115746_complete_metadata_28706-cells.csv.gz")

if (!file.exists(exon_path)) {
  system(sprintf('curl --retry 10 --retry-delay 5 --ftp-method nocwd -o %s "%s"', exon_path, exon_url))
  system(sprintf("gzip -t %s", exon_path))  # MUST pass -- see download lesson above
}
if (!file.exists(meta_path)) {
  system(sprintf('curl -L --retry 5 --retry-delay 5 -o %s "%s"', meta_path, meta_url))
}

meta <- fread(meta_path)
real_cells <- meta[organism == "Mus musculus" & !(cell_class %in% c("Low Quality", "No Class"))]
cat("Filtered real cells in metadata:", nrow(real_cells), "\n")

# Chunked read + sparsify -- see memory lesson above. Do not read in one pass.
header <- fread(exon_path, nrows = 0)
cell_ids <- colnames(header)[-1]
n_genes <- 45768
chunk_size <- 5000
chunks <- list()
gene_ids_all <- character(0)
start <- 1; i <- 1
while (start <= n_genes) {
  end <- min(start + chunk_size - 1, n_genes)
  chunk <- fread(exon_path, skip = start, nrows = end - start + 1, header = FALSE)
  gene_ids_all <- c(gene_ids_all, chunk[[1]])
  mat_chunk <- as.matrix(chunk[, -1, with = FALSE])
  rm(chunk); gc(verbose = FALSE)
  chunks[[i]] <- as(mat_chunk, "CsparseMatrix")
  rm(mat_chunk); gc(verbose = FALSE)
  start <- end + 1; i <- i + 1
}
counts_mat <- do.call(rbind, chunks)
rm(chunks); gc(verbose = FALSE)
rownames(counts_mat) <- gene_ids_all
colnames(counts_mat) <- cell_ids

common_ids <- intersect(colnames(counts_mat), real_cells$sample_name)
cat("Matched cells (metadata x counts):", length(common_ids), "\n")

counts_filtered <- counts_mat[, common_ids]
rm(counts_mat); gc(verbose = FALSE)
meta_filtered <- real_cells[match(common_ids, real_cells$sample_name)]
stopifnot(identical(meta_filtered$sample_name, colnames(counts_filtered)))

sce <- SingleCellExperiment(
  assays = list(counts = counts_filtered),
  colData = DataFrame(
    sample_name = meta_filtered$sample_name,
    cell_class = meta_filtered$cell_class,
    cell_subclass = meta_filtered$cell_subclass,
    cell_cluster = meta_filtered$cell_cluster,
    dissected_region = meta_filtered$dissected_region,
    donor_id = meta_filtered$donor_id,
    donor_sex = meta_filtered$donor_sex
  )
)

out_path <- "/mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_sce.rds"
saveRDS(sce, out_path, compress = TRUE)
cat("Saved:", out_path, "\n")
cat("Final dimensions:", dim(sce), "\n")
