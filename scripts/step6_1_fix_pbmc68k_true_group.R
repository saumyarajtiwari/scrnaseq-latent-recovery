# Fix: propagate the corrected PBMC68k re-annotation labels into the 6
# pbmc68k embedding .rds files, whose internally-stored true_group field
# was never updated when embedding_manifest.csv was re-joined (metadata-only)
# after the original re-annotation fix (see project history: "embedding_manifest.csv
# re-joined (not re-scanned) to propagate the fix cheaply, since only metadata
# changed, not file integrity" — that re-join never touched the .rds files themselves).
#
# ALIGNMENT CAVEAT (logged into each file below): no independent per-cell ID
# exists in either the source SCE (no colnames) or the embedding object to
# directly verify row order. This fix assumes positional alignment holds,
# based on: (1) cell counts match exactly (65,690 = 65,690), and (2) the
# project history describes the re-annotation as relabeling cells on the
# existing harmonized file in place, not reordering/resubsetting them.
# This is a documented assumption, not a proven alignment (unlike the Baron
# check earlier in Step 6.1, which had an independent batch_id cross-check).

suppressPackageStartupMessages(library(SingleCellExperiment))

src_path <- "/mnt/extra/scrnaseq-data/real/pbmc68k/pbmc68k_harmonized.rds"
src <- readRDS(src_path)
src_true_group <- as.character(colData(src)$true_group)

cat("Source true_group distribution:\n")
print(table(src_true_group, useNA = "always"))

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
pbmc_files <- manifest[manifest$source == "pbmc68k", ]
cat(sprintf("\nFound %d pbmc68k embedding files to check/fix\n", nrow(pbmc_files)))

alignment_note <- paste0(
  "true_group was stale (100% NA) from before PBMC68k re-annotation; patched ",
  Sys.Date(), " from corrected source SCE (5 canonical types, 65690/65690 labeled). ",
  "ASSUMPTION: positional cell-order alignment between source and embedding ",
  "could not be independently verified (no cell_id/colnames in either object); ",
  "based on matching cell counts and documented in-place relabeling (no ",
  "reorder/resubset) per project history. Embedding/loadings unchanged."
)

for (i in seq_len(nrow(pbmc_files))) {
  fp <- pbmc_files$file_path[i]
  cat(sprintf("\n--- %s (%s) ---\n", pbmc_files$method[i], fp))

  obj <- readRDS(fp)

  n_obj <- length(obj$true_group)
  if (n_obj != length(src_true_group)) {
    stop(sprintf("ABORT: cell count mismatch (%d in embedding vs %d in source) for %s",
                  n_obj, length(src_true_group), fp))
  }

  old_na_count <- sum(is.na(obj$true_group))
  cat(sprintf("  Before fix: %d/%d NA in true_group (class: %s)\n",
              old_na_count, n_obj, class(obj$true_group)))

  # Preserve original type (factor) explicitly rather than relying on
  # implicit coercion
  obj$true_group <- factor(src_true_group)

  new_na_count <- sum(is.na(obj$true_group))
  cat(sprintf("  After fix:  %d/%d NA in true_group (class: %s)\n",
              new_na_count, n_obj, class(obj$true_group)))
  print(table(obj$true_group, useNA = "always"))

  if (is.null(obj$flags)) obj$flags <- list()
  obj$flags$true_group_patched <- alignment_note

  tmp_path <- paste0(fp, ".tmp")
  saveRDS(obj, tmp_path)
  file.rename(tmp_path, fp)
  cat("  Written (atomic).\n")
}

cat("\n=== Re-verification of all 6 files ===\n")
for (i in seq_len(nrow(pbmc_files))) {
  obj <- readRDS(pbmc_files$file_path[i])
  cat(sprintf("%s -> n_NA: %d / %d | class: %s | flag present: %s\n",
              pbmc_files$method[i], sum(is.na(obj$true_group)), length(obj$true_group),
              class(obj$true_group), !is.null(obj$flags$true_group_patched)))
}
