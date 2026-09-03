# EDA Checkpoint 1 — Code Gap Documentation

## What happened

EDA Checkpoint 1 (simulated-output verification: sparsity, depth, and
DE-magnitude checks, plus a summary pass/fail table) was carried out as a
series of interactive, ad hoc R commands during a working session, rather
than as a saved, standalone `.R` script. Its outputs — the tables and figures
listed below — were committed to this repository; the exact code that
produced them was not.

This is a genuine gap, not an oversight being hidden: it is called out
explicitly here, in `PATCH_NOTES.md`, and in `PROJECT_HANDOVER.md`, and the
tracked outputs below can be independently sanity-checked against the
simulation grid at any time.

## What is tracked and available right now

- `results/step1_eda/tables/eda_1_1_*.csv` through `eda_1_4_pass_fail_table.csv`
- `results/step1_eda/figures/eda_1_*.png`
- A narrative summary of the findings in `PROJECT_HANDOVER.md`, Step 1.7/1.8 sections.

## What the checkpoint did, in enough detail to reconstruct it

Per the project record, EDA Checkpoint 1 performed four checks against the
already-generated and already-validated simulation grid
(`data/simulated/{splatter,scdesign3,symsim}/`, using
`data/simulated/param_grid.csv` for expected values):

1. **Sparsity verification (1.1)**: recomputed actual sparsity per file
   directly from raw counts (`1 - nnzero(counts) / length(counts)`) and
   compared against each file's stored `run_params$sparsity_label` /
   `actual_sparsity`, following the same recomputation logic later reused
   and extended in `validate_output_inventory.R` (Step 1.7).
2. **Depth verification (1.2)**: recomputed mean per-cell UMI depth
   (`mean(colSums(counts))`) and compared against each simulator's own
   calibration table (`{simulator}_calib_depth.csv`), the same comparison
   later formalized in Step 1.7's `get_expected_depth()` function.
3. **DE-magnitude check (1.3)**: computed log2 fold-change between true
   biological groups using each file's `true_group_means` (or an equivalent
   early version of it), following the identical `log2(mean+1)` formula
   later reused directly in Step 6.5b (`step6_5b_simulated_log2fc_discovery.R`,
   which explicitly notes it reuses "the same formula as the original EDA
   1.3 investigation").
4. **Summary pass/fail rollup (1.4)**: aggregated the above three checks per
   simulator/parameter-combination into a single pass/fail table, following
   the same per-group aggregation pattern used in Step 1.7's
   `validation_summary.csv`.

## Recommended path to closing this gap

The cleanest reconstruction is to **adapt `step6_5b_simulated_log2fc_discovery.R`'s
DE-magnitude logic and `validate_output_inventory.R`'s sparsity/depth
recomputation logic** into a single new script
(`code/eda_checkpoint1/eda_checkpoint1_reconstructed.R`), since both already
contain verified, working implementations of the exact same computations
this checkpoint originally performed. This has not yet been done. Until it
is, this checkpoint should be described in the manuscript as "verification
outputs are available; the exact reconstructing script is provided as
`eda_checkpoint1_reconstructed.R`, built from the equivalent, already-verified
logic used in later pipeline stages" once written — or, if left unreconstructed,
described plainly as "performed via interactive analysis; outputs retained,
exact commands not preserved as a standalone script," which is the accurate
statement of the current state.
