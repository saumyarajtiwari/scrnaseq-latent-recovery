# Real Dataset Metadata Catalog

Step 2.7. This document is the authoritative reference for cell-type labels,
batch structure, tissue origin, sequencing technology, and known biological
caveats for all six real datasets used in this project. All downstream
latent-recovery metrics computed on real data should be interpreted against
the ground truth and limitations documented here, not assumed.

All cell/gene/label counts below are extracted directly from the final
harmonized SCE objects (Step 2.6 output), not transcribed from memory.

---

## 1. PBMC 68k

- **Source**: 10x Genomics "Fresh 68k PBMCs (Donor A)", via `TENxPBMCData`
  Bioconductor package (Step 2.1).
- **Tissue of origin**: Human peripheral blood mononuclear cells (PBMC),
  single donor ("Donor A").
- **Sequencing technology**: 10x Genomics droplet-based 3' scRNA-seq.
- **Post-QC dimensions**: 20,387 genes x 65,690 cells.
- **Cell-type labels**: **NONE AVAILABLE.** This dataset ships with no
  manually-curated or bundled cell-type annotation. The only official
  "annotation" elsewhere for this dataset is a classifier-inferred soft
  label (correlation against 10 FACS-purified reference populations) plus
  an incomplete/inconsistent k-means cluster file (40k of 68k cells).
  Decision (Step 2.1): do not use either; re-annotate independently later
  via a marker-based approach consistent with the PBMC 3k method (Step 1.4).
  **This re-annotation has not yet been performed as of Step 2.7.**
  `true_group` is NA for all 65,690 cells in the harmonized object.
- **Batch structure**: None (single donor, single sample; `n_batches = 1`,
  trivial).
- **Known biological structure**: Not yet assigned. Literature expectation
  for PBMC is ~8-10 major immune populations (T cell subsets, B cells, NK
  cells, monocytes, dendritic cells) but this should not be treated as
  ground truth until actual re-annotation is performed.
- **Sparsity (post-QC)**: 97.3%.

---

## 2. Muraro (GSE85241)

- **Source**: `scRNAseq::MuraroPancreasData()` (Step 2.2).
- **Tissue of origin**: Human pancreas.
- **Sequencing technology**: CEL-seq2 (plate-based).
- **Post-QC dimensions**: 18,197 genes x 2,403 cells.
- **Batch structure**: 4 donors (D28, D29, D30, D31) — original plate-based
  batch identity, preserved unchanged.
- **Cell-type labels**: 10 published types (acinar, alpha, beta, delta,
  duct, endothelial, epsilon, mesenchymal, pp) plus 277 cells (11.5%) with
  no author-assigned label (`true_group = NA`).
- **Known caveat**: QC (Step 2.5) confirmed that ALL cells removed during
  quality filtering came from the unlabeled category — 0% of any named
  cell type was affected by QC. The 277 remaining unlabeled cells are
  QC-passing but genuinely un-annotated by the original authors, not a
  processing artifact on our end.
- **Sparsity (post-QC)**: 72.0%.

---

## 3. Baron (GSE84133, human)

- **Source**: `scRNAseq::BaronPancreasData(which = "human")` (Step 2.2).
- **Tissue of origin**: Human pancreas.
- **Sequencing technology**: inDrop (droplet-based).
- **Post-QC dimensions**: 17,499 genes x 8,569 cells.
- **Batch structure**: 4 donors (GSM2230757-760).
- **Cell-type labels**: 14 published types, fully labeled (0 unlabeled
  cells) — acinar, activated_stellate, alpha, beta, delta, ductal,
  endothelial, epsilon, gamma, macrophage, mast, quiescent_stellate,
  schwann, t_cell. Includes rare populations (schwann n=13, t_cell n=7 in
  raw data) — QC confirmed 0% cells flagged in this dataset overall, so
  these rare populations are fully intact post-QC.
- **Sparsity (post-QC)**: 89.2%.

---

## 4. Segerstolpe (E-MTAB-5061)

- **Source**: `scRNAseq::SegerstolpePancreasData()` (Step 2.2).
- **Tissue of origin**: Human pancreas.
- **Sequencing technology**: Smart-seq2 (full-length, plate-based).
- **Post-QC dimensions**: 23,238 genes x 2,929 cells.
- **Batch structure**: 10 individuals. One (`AZ`) has only 96 cells
  pre-QC vs. 383-384 for the other nine — a real, published asymmetry,
  not a download defect.
- **Disease covariate**: 10 individuals include both normal and Type II
  diabetes mellitus donors, confounded with donor identity (this is real
  biological structure in the source data, not introduced during
  processing). Preserved as an explicit extra `colData` field (`disease`).
- **Cell-type labels**: 14 published types post-QC, plus 720 cells (24.6%)
  with no author-assigned label. Note: original metadata also contained
  small explicit "unclassified cell" (n=2) and "unclassified endocrine
  cell" (n=41) categories, distinct from the larger blank/NA group.
- **Known caveat**: Same pattern as Muraro — QC (Step 2.5) confirmed 100%
  of QC-removed cells came from the unlabeled group; 0% impact on any
  named cell type.
- **Sparsity (post-QC)**: 77.1%.

---

## 5. Tabula Sapiens — Lung subset

- **Source**: CZ CELLxGENE Census, `cellxgene.census` R package, pinned to
  Census version `2025-11-08` (Step 2.3). `dataset_id =
  0d2ee4ac-05ee-40b2-afb6-ebb584caa867`, confirmed as the unique Tabula
  Sapiens Lung entry before any cell data was pulled.
- **Tissue of origin**: Human lung.
- **Sequencing technology**: **Mixed within this single dataset** — both
  Smart-seq2 and 10x 3' v3. This constitutes an internal technology batch
  axis distinct from donor identity, preserved as an explicit extra
  `colData` field (`assay`).
- **Post-QC dimensions**: 56,139 genes x 61,292 cells.
- **Batch structure**: 4 donors, plus the above technology-based sub-batch
  structure.
- **Cell-type labels**: 34 published types, fully labeled (0 unlabeled
  cells). Deliberately selected (Step 2.3, per original project scope) for
  its broad, sparse spread of many distinct types rather than a small
  number of dense, well-separated clusters.
- **Sparsity (post-QC)**: 93.0%.

---

## 6. Tasic et al. 2018 — Mouse Cortex (GSE115746)

- **Source**: Direct GEO/FTP download (no Bioconductor wrapper exists for
  this specific dataset; NOT to be confused with `scRNAseq::TasicBrainData()`,
  which is the earlier, unrelated Tasic et al. 2016 dataset) (Step 2.4).
- **Tissue of origin**: Mouse visual cortex (VISp) and anterior lateral
  motor cortex (ALM) — preserved as extra `colData` field
  (`dissected_region`).
- **Sequencing technology**: Smart-seq2 (full-length, plate-based).
- **Post-QC dimensions**: 42,865 genes x 20,562 cells.
- **Batch structure**: **341 unique donor IDs** (individual mice, confirmed
  as real Allen Institute animal identifiers spanning many transgenic
  Cre-driver lines; median 32 cells/animal, range 1-334). This is a
  fundamentally different, much finer-grained batch-complexity profile
  than every other real dataset here (4-10 batches each) — reflects this
  study's actual experimental design and should be treated as such, not
  normalized or treated as anomalous.
- **Cell-type label hierarchy** (three levels, all preserved):
  - `true_group` / **top-level, used for primary analysis**: `cell_subclass`,
    24 categories. **Deviates from original project scope**, which
    referenced an unconfirmed "11 top-level classes" scheme that could not
    be traced to anything in the actual published metadata (neither the
    4-category `cell_class` nor the 24-category `cell_subclass` equals 11).
    Decision (Step 2.4, delegated project judgment): use `cell_subclass`
    as top-level, since it is an official, citable Allen Institute
    metadata column.
  - `cell_class`: 4 categories (Endothelial, GABAergic, Glutamatergic,
    Non-Neuronal) — used for QC blocking (Step 2.5), NOT as the primary
    top-level label.
  - `cell_cluster`: 134 categories (fine-grained; original scope referenced
    133 — off by one, not individually root-caused, consistent with an
    unclassified catch-all appearing elsewhere in this same metadata).
- **Known caveats**:
  - **Cell-count discrepancy**: original scope referenced 23,822 cells;
    actual post-QC count is 20,562. Root cause chain (Step 2.4): raw
    metadata had 28,706 total sample records (2,986 controls/spike-ins
    excluded) -> 24,061 real, non-excluded mouse cells -> only 22,113 of
    those had matching columns in the exon-counts matrix (1,948-cell gap,
    investigated and found to be proportionally spread across nearly every
    subclass and both regions, all QC-passing — most plausibly a data-
    availability gap in this specific GEO release, not a targeted
    exclusion) -> 20,562 after QC (Step 2.5) removed a further 1,551 cells
    via cell_class-blocked adaptive thresholds.
  - **CR (Cajal-Retzius) population**: independently flagged as a likely
    non-cortical contaminant population in Step 2.4, and this is
    independently supported by Step 2.5's QC results, where CR showed the
    highest flagged fraction (88.2%) even after correct blocking — two
    independent lines of evidence pointing the same direction.
  - **Meis2 population**: small (n=55 pre-QC), showed elevated flagged
    fraction (74.5%) post-blocking; likely reflects genuine within-group
    variance at low sample size rather than a blocking failure, but not
    independently confirmed further.
- **Sparsity (post-QC)**: 78.2%.

---

## Summary Table

| Dataset | Tissue | Technology | Cells (post-QC) | Genes (post-QC) | Cell types | Batches | Fully labeled? |
|---|---|---|---|---|---|---|---|
| pbmc68k | Human PBMC | 10x 3' droplet | 65,690 | 20,387 | 0 (deferred) | 1 | No — 100% unlabeled |
| muraro | Human pancreas | CEL-seq2 | 2,403 | 18,197 | 10 | 4 | No — 11.5% unlabeled |
| baron | Human pancreas | inDrop | 8,569 | 17,499 | 14 | 4 | Yes |
| segerstolpe | Human pancreas | Smart-seq2 | 2,929 | 23,238 | 14 | 10 | No — 24.6% unlabeled |
| ts_lung | Human lung | Smart-seq2 + 10x 3'v3 | 61,292 | 56,139 | 34 | 4 | Yes |
| tasic2018 | Mouse cortex | Smart-seq2 | 20,562 | 42,865 | 24 | 341 | Yes |
