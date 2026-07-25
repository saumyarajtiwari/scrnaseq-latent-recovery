# EDA 2.3 — Batch Effect Visibility Check (Pancreas Only): Findings

## What was run
Three PCA analyses on the three pancreas datasets (Muraro, Baron human,
Segerstolpe), gene set restricted to 15,511 genes common to all three
(Muraro's `SYMBOL__chrN` rowname suffix stripped first to enable matching):

1. **Combined, raw counts** (as literally specified by this step) — `prcomp()`
   on the full dense matrix (13,901 cells x 15,511 genes).
2. **Combined, log1p counts** (supplementary, added after #1 proved
   uninformative) — `irlba` truncated SVD on the sparse matrix directly.
3. **Per-dataset, log1p counts, donor-colored** (closing check, added after
   #2 also proved uninformative for the step's actual question) — `irlba`
   run separately within each of the three datasets.

## Results

**#1 (combined, raw):** PC1 = 57.4%, PC2 = 14.9%. Segerstolpe's mean
per-cell depth (~375,000 UMIs, ~64x higher than Baron's ~5,800) completely
dominated the plot scale; Muraro and Baron were visually compressed into a
single indistinguishable point at the origin. Not informative for
assessing batch or biological structure in the two lower-depth datasets.

**#2 (combined, log1p):** PC1 = 82%, PC2 = 7% — MORE concentrated in PC1
than the raw-count version, not less. Produced a classic PCA "horseshoe"
(Guttman) artifact, a known signature of a single dominant gradient (here,
still depth/technology, since Segerstolpe alone spans a large internal
depth range — recall its `AZ` donor has only 96 cells vs. 383-384 for the
other nine). Donor colors showed some blocky segments along the arc for
Segerstolpe specifically, but this is not distinguishable from donors
simply having different characteristic depths and landing at different
points along the same depth-driven curve — not confirmed as genuine
biological/batch clustering.

**#3 (per-dataset, log1p, donor-colored):** Donor colors were broadly
INTERMIXED in all three datasets (Muraro: PC1=60.1%/PC2=17.9%; Baron:
PC1=38%/PC2=24.4%; Segerstolpe: PC1=70.1%/PC2=13.6%). No dataset showed
donors forming visually separated clusters. Each again showed an
arc/horseshoe-like structure, likely reflecting cell-type identity as the
dominant gradient rather than donor batch identity.

## Interpretation
This step's stated goal was to confirm that batch structure is visually
detectable in these three pancreas datasets, validating their use for
batch-complexity stress testing in Step 3. **This could not be confirmed
via simple raw- or log1p-count PCA, at either the combined or per-dataset
level.** This is not evidence that batch effects are absent — donor-level
technical variation in this exact three-dataset pancreas collection is
well-established in the integration-benchmarking literature — but it does
mean naive PCA on unprocessed counts is not a sufficient diagnostic to
demonstrate it. Cell-type-driven biological variation appears to dominate
the first few principal components in all three datasets, likely masking
weaker batch-driven variation at this level of analysis.

## Scope implication (for manuscript)
Do not claim, based on this EDA checkpoint alone, that batch structure was
confirmed as visually detectable in raw form. The more accurate framing is:
raw/log1p PCA was tested and found insufficient to reveal batch structure,
motivating (rather than validating in advance) the need for the
normalization, HVG selection, and/or explicit batch-integration methods
planned in Step 3 before batch-complexity stress-testing can be properly
evaluated. This should be revisited after Step 3's normalization pipeline
is in place, ideally with a dedicated re-check of batch visibility on
processed data.
