#!/bin/bash
# Step 2.6 — Format Harmonization driver for all six real datasets.
#
# Target schema matched to Step 1's simulated SCEs where a real equivalent
# exists (assay="counts" only, raw/unnormalized -- Step 3 handles
# normalization). Simulated-only fields (run_id, true_group_means,
# ground_truth_source="simulator", simulator name) have NO defensible real-
# data equivalent and are NOT force-mapped:
#   - true_group_means would require clustering/embedding the data first,
#     which is circular (that's the method under evaluation, not a ground
#     truth input).
#   - run_id/simulator are simulation-only concepts.
#
# Real-dataset-specific colData mapping used:
#   pbmc68k:      true_group=NA (no annotation exists; Step 2.1 decision,
#                 re-annotation deferred), batch_id=Sample
#   muraro:       true_group=label, batch_id=donor
#   baron:        true_group=label, batch_id=donor
#   segerstolpe:  true_group=`cell type`, batch_id=individual,
#                 extra: disease (T2D covariate, confounded with donor --
#                 real biology, not introduced here)
#   ts_lung:      true_group=cell_type, batch_id=donor_id,
#                 extra: assay (Smart-seq2 vs 10x 3'v3 -- internal technology
#                 batch axis within this single dataset)
#   tasic2018:    true_group=cell_subclass (24 cats, per Step 2.4's
#                 documented decision to use this over the unconfirmed
#                 11-class spec), batch_id=donor_id,
#                 extra: cell_class, cell_cluster, dissected_region
#
# BUG FOUND AND FIXED: pbmc68k's SCE (from TENxPBMCData) has EMPTY
# colnames() -- length 0, not ncol(sce). The real per-cell barcode lives in
# the `Barcode` colData column instead. harmonize_one_dataset.R now checks
# colnames() length and falls back to `Barcode` when empty, with an explicit
# stopifnot() length check rather than silently producing a malformed
# DataFrame (which is what happened on first attempt -- "different row
# counts implied by arguments" error).
#
# NOTABLE FINDING (not a bug): tasic2018's donor_id has 341 unique values
# (median 32 cells/donor, range 1-334) -- confirmed as real Allen Institute
# animal IDs, consistent with this study's design (many transgenic mouse
# lines, few cells per animal). This is a genuinely different batch-
# complexity profile than the other five real datasets (4-10 batches each)
# and should be treated as real experimental structure, not normalized away
# or treated as an anomaly.
#
# Each dataset processed as a SEPARATE Rscript process (memory safety,
# per convention since Step 2.3). Raw and QC-stage intermediate files were
# deleted after harmonized versions were verified (both stages fully
# reproducible via download_*.R / qc_filter_one_dataset.R scripts).

set -e
cd ~/Desktop/scrnaseq-latent-recovery

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/pbmc68k/pbmc68k_qc.rds \
  /mnt/extra/scrnaseq-data/real/pbmc68k/pbmc68k_harmonized.rds \
  pbmc68k NA Sample NA

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/muraro/muraro_qc.rds \
  /mnt/extra/scrnaseq-data/real/muraro/muraro_harmonized.rds \
  muraro label donor NA

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/baron/baron_qc.rds \
  /mnt/extra/scrnaseq-data/real/baron/baron_harmonized.rds \
  baron label donor NA

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_qc.rds \
  /mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_harmonized.rds \
  segerstolpe "cell type" individual disease

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/tabula_sapiens_lung/tabula_sapiens_lung_qc.rds \
  /mnt/extra/scrnaseq-data/real/tabula_sapiens_lung/tabula_sapiens_lung_harmonized.rds \
  ts_lung cell_type donor_id assay

Rscript R/02_real_data/harmonize_one_dataset.R \
  /mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_qc.rds \
  /mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_harmonized.rds \
  tasic2018 cell_subclass donor_id "cell_class,cell_cluster,dissected_region"

echo "NOTE: raw and QC-stage intermediate files are deleted after this script"
echo "runs successfully and outputs are verified -- see Step 2.6 commit for"
echo "the deletion commands used. Re-run download_*.R + qc_filter_one_dataset.R"
echo "first if intermediates are needed again."
