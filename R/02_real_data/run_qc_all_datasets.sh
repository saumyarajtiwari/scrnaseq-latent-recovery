#!/bin/bash
# Step 2.5 — QC filtering driver for all six real datasets.
# Each dataset invoked as a SEPARATE Rscript process (memory safety, per
# project convention established in Steps 2.3-2.4).
#
# QC approach: adaptive MAD-based thresholds (scuttle::isOutlier, 3 MADs,
# log-scale for sum/detected, one-sided) rather than fixed universal cutoffs,
# per this step's explicit requirement not to over-filter real rare
# populations.
#
# Mitochondrial-gene availability investigated per dataset (see commit
# history / handover for full diagnostic trail): only pbmc68k (via
# ENSEMBL_ID) and tabula_sapiens_lung (via Ensembl rowname) have detectable
# canonical mito genes in their gene panels. muraro, baron, segerstolpe, and
# tasic2018 genuinely lack mitochondrial genes in their published gene sets
# (upstream exclusion by original authors) -- mito% is not used as a QC
# covariate for these four; this is a documented dataset limitation, not an
# oversight.
#
# CRITICAL FINDING: global (unblocked) MAD thresholds on tasic2018 would
# remove 90-100% of every non-neuronal cell type (Peri, Macrophage, SMC,
# Astro, Endo, VLMC) due to genuinely lower RNA content in glia/vasculature
# vs. neurons under Smart-seq2 -- NOT a quality problem. Fixed by blocking
# isOutlier() by cell_class (4 groups), the standard OSCA-recommended
# approach for datasets with known, biologically distinct cell populations.
# This reduced tasic2018's flagged fraction from 11.7% (global, biased) to
# 7.0% (blocked, biologically appropriate) while protecting every
# non-neuronal population. Remaining higher-than-average flagged fractions
# in CR (88.2%) and Meis2 (74.5%) after blocking are consistent with CR's
# already-flagged likely-contaminant status (Step 2.4) and Meis2's small
# sample size (n=55), not a blocking failure.
#
# muraro and segerstolpe: confirmed via per-cell-type breakdown that 100%
# of flagged cells fall into the authors' own unclassified/blank-label
# category -- 0% of any named, annotated cell type is affected. No blocking
# needed for these two.

set -e
cd ~/Desktop/scrnaseq-latent-recovery

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/pbmc68k/pbmc68k_raw.rds \
  /mnt/extra/scrnaseq-data/real/pbmc68k/pbmc68k_qc.rds \
  ensembl_pbmc68k NA

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/muraro/muraro_raw.rds \
  /mnt/extra/scrnaseq-data/real/muraro/muraro_qc.rds \
  none NA

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/baron/baron_raw.rds \
  /mnt/extra/scrnaseq-data/real/baron/baron_qc.rds \
  none NA

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/segerstolpe/segerstolpe_raw.rds \
  /mnt/extra/scrnaseq-data/real/segerstolpe/segerstolpe_qc.rds \
  none NA

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/tabula_sapiens_lung/tabula_sapiens_lung_raw.rds \
  /mnt/extra/scrnaseq-data/real/tabula_sapiens_lung/tabula_sapiens_lung_qc.rds \
  ensembl_tslung NA

Rscript R/02_real_data/qc_filter_one_dataset.R \
  data/real/tasic2018/tasic2018_sce.rds \
  /mnt/extra/scrnaseq-data/real/tasic2018/tasic2018_qc.rds \
  none cell_class
