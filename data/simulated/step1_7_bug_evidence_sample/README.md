# Step 1.7 Bug Evidence Sample

Small representative sample (3 fit-groups x 5 sparsity labels = 15 files
per simulator) extracted from the full original pre-fix scDesign3/SymSim
backups before the 28GB bulk backups were deleted, per Step 1.8 disk-space
needs.

Each group holds every real parameter fixed (depth, dropout, separability,
n_cells, batch, gene_strategy, clipping) and varies only sparsity_label --
directly reproducing Step 1.7's original finding that sparsity_label had
zero generative effect (byte-identical matrices across all 5 labels within
a fit-key group). See PROJECT_HANDOVER.md for the full Step 1.7 narrative.

Full bulk backups (28GB) deleted after this sample was extracted and
verified to still reproduce the original finding.
