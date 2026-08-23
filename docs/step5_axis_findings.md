# Step 5.2–5.7: Axis-Level Findings and Limitations

Documented during visual review of the six single-axis subspace-recovery
plots. Each finding below was investigated with direct data evidence,
not assumption, following this project's established diagnostic standard.

## 1. Dropout — scDesign3 label swap (FIXED)

See `scripts/step5_fix_scdesign3_dropout_labels.R` for full detail.
scDesign3's "low"/"high" dropout labels were found swapped relative to
their intended severity (confirmed via achieved_sparsity vs. label
across 6 parameter contexts, 6/6 consistent). Fixed at the results-table
level (not param_grid.csv, which is simulator-shared and correctly
specifies the intended design). SymSim's separate, non-swap-fixable
dropout non-monotonicity remains unresolved (see #3 below).

## 2. Clipping — confirmed generatively inert for scDesign3 and SymSim (NOT FIXED - documented limitation)

Visual review of Step 5.7's plot showed scDesign3's and SymSim's method
curves completely flat across all three clipping levels (none,
log_stabilized, clip99), while Splatter showed genuine variation.

Confirmed via direct inspection of raw count matrices at matched
parameter contexts: scDesign3's max_count and total_counts were
byte-for-byte IDENTICAL (367, 720755) across all three clipping labels.
SymSim showed the same pattern (899, 1372987, identical across all
three). Splatter showed real differences (320/1138283 vs 183/1149711
vs 247/1141012), confirming its clipping mechanism functions correctly.

This is the same class of defect already documented in Step 1.7 for
scDesign3's and SymSim's sparsity/dropout axes ("generatively inert"),
now confirmed as a third instance affecting the clipping axis for both
simulators. Unlike the dropout finding, this is not a simple label
swap and has no cheap metadata-only fix - the underlying simulated
data genuinely does not respond to the clipping parameter for these
two simulators.

**Manuscript treatment: Splatter should be treated as the sole
reliable simulator for clipping-axis conclusions.** scDesign3/SymSim's
flat clipping curves should be reported as a confirmed simulator
limitation, not as evidence that clipping has no effect on recovery
in general.

## 3. Batch complexity — SymSim shows an unexplained positive (backwards) relationship (NOT FIXED - documented limitation)

Visual review of Step 5.5's plot showed SymSim's Raw PCA and Libnorm
PCA recovery increasing sharply (~0.44 to ~0.73) as batch complexity
increased from none to complex - the opposite of the expected
direction (more batch complexity should generally make recovery
harder, not easier). Splatter and scDesign3 did not show this pattern.

Two standard explanations were tested directly and both ruled out:

- **Batch/group label confounding**: cross-tabulated batch_id against
  true_group for SymSim's complex-batch run. Perfectly balanced
  (100/100/100/100/100 in batch 1, 60/60/60/60/60 in batch 2, etc.) -
  no correlation between batch assignment and biological group.
  Splatter showed the same balanced pattern. Confounding ruled out.
- **Batch-effect inertness** (i.e., the same class of bug as clipping
  above): compared mean library size and per-gene mean expression
  across batches directly. SymSim showed real, substantial per-gene
  shifts between batches (e.g., one gene's mean expression varied
  0.64 -> 6.17 -> 0.14 across the three batches, roughly a 40x swing) -
  larger in magnitude than Splatter's corresponding shifts on the same
  genes. The batch effect is confirmed real and substantial, not inert.

With both standard explanations ruled out via direct evidence, the
remaining candidate explanation is a more subtle interaction between
SymSim's specific batch-effect generation mechanism and its biological
signal structure - understanding this precisely would require reading
and analyzing SymSim's internal batch-effect generation source code,
a materially different and deeper investigation than the data-level
diagnostics used elsewhere in this step, and judged disproportionate
to Step 5's scope (downstream analysis, not simulator development).

**Manuscript treatment: report as an open, confirmed-real, ruled-out-
the-obvious-causes finding.** State plainly that confounding and
inertness were tested and excluded, and that the precise generative
mechanism is left to future investigation, consistent with standard
practice for documenting unexplained-but-verified anomalies.
