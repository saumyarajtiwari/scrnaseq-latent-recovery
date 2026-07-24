# Step 2.2 — Pancreas Datasets Download
# Downloads Muraro (GSE85241), Baron (GSE84133, human only), and
# Segerstolpe (E-MTAB-5061) via the Bioconductor scRNAseq package,
# each backed by ExperimentHub with local caching.
#
# These three datasets deliberately differ in sequencing technology
# (CEL-seq2, inDrop, Smart-seq2 respectively) and are used for
# batch-complexity stress testing. Original donor/individual batch
# identity labels and published cell-type annotations are preserved
# unchanged from each dataset's own colData.
#
# Segerstolpe additionally carries a disease covariate (normal vs.
# type II diabetes mellitus) confounded with donor identity — this
# is real biological structure in the source data, not introduced
# here, and should be treated as a known covariate in any batch
# analysis using this dataset.

library(scRNAseq)
library(SingleCellExperiment)
library(Matrix)

save_verified <- function(sce, out_path, batch_col, label_col) {
  cts <- assay(sce, "counts")
  if (!is(cts, "dgCMatrix")) {
    counts(sce) <- as(cts, "dgCMatrix")
  }
  saveRDS(sce, out_path, compress = TRUE)

  check <- readRDS(out_path)
  ok <- identical(dim(check), dim(sce)) &&
        is(assay(check, "counts"), "dgCMatrix") &&
        identical(table(check[[batch_col]]), table(sce[[batch_col]])) &&
        identical(table(check[[label_col]]), table(sce[[label_col]]))

  cat(out_path, "- verified:", ok, "\n")
  if (!ok) stop("Verification failed for ", out_path)
}

# --- Muraro (GSE85241, CEL-seq2) ---
muraro <- MuraroPancreasData()
save_verified(
  muraro,
  "~/Desktop/scrnaseq-latent-recovery/data/real/muraro/muraro_raw.rds",
  batch_col = "donor", label_col = "label"
)

# --- Baron (GSE84133, inDrop, human only) ---
baron <- BaronPancreasData(which = "human")
save_verified(
  baron,
  "~/Desktop/scrnaseq-latent-recovery/data/real/baron/baron_raw.rds",
  batch_col = "donor", label_col = "label"
)

# --- Segerstolpe (E-MTAB-5061, Smart-seq2) ---
seger <- SegerstolpePancreasData()
save_verified(
  seger,
  "~/Desktop/scrnaseq-latent-recovery/data/real/segerstolpe/segerstolpe_raw.rds",
  batch_col = "individual", label_col = "cell type"
)

cat("\nAll three pancreas datasets downloaded, saved, and verified.\n")
