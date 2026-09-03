#!/bin/bash
set -x

echo "=== [1/8] Step 6.3 discovery: $(date) ==="
Rscript scripts/step6_3_discover_null_control_files.R

echo "=== [2/8] Step 6.3 phantom clustering: $(date) ==="
Rscript scripts/step6_3_phantom_clustering_v2.R

echo "=== [3/8] Step 6.4 variance hijacking: $(date) ==="
Rscript scripts/step6_4_variance_hijacking.R

echo "=== [4/8] Step 6.5a pancreas DPT: $(date) ==="
Rscript scripts/step6_5a_pancreas_dpt.R

echo "=== [5/8] Step 6.5b log2FC discovery: $(date) ==="
Rscript scripts/step6_5b_simulated_log2fc_discovery.R

echo "=== [6/8] Step 6.5c variance ratio: $(date) ==="
Rscript scripts/step6_5c_simulated_variance_ratio.R

echo "=== [7/8] Step 6.6 prebuild + neighborhood collapse: $(date) ==="
Rscript scripts/step6_6_prebuild_real_reference_spaces.R
Rscript scripts/step6_6_neighborhood_collapse.R

echo "=== [8/8] Step 6.7, 6.7b, 6.8: $(date) ==="
Rscript scripts/step6_7_subspace_rotation_slippage.R
Rscript scripts/step6_7b_slippage_gradient_analysis.R
Rscript scripts/step6_8_cross_reference_boundaries.R

echo "=== ALL STEPS COMPLETE: $(date) ==="
