#!/bin/bash
# =============================================================================
# migrate_repo_structure.sh
#
# Reorganizes scrnaseq-latent-recovery into a clean, publication-ready
# structure, entirely for GitHub (no Zenodo/external hosting steps here).
#
# SAFETY DESIGN:
#   - Refuses to run unless `git status` is clean (no uncommitted changes)
#   - Creates a backup branch before touching anything
#   - Uses `git mv` throughout, so full history is preserved on every file
#   - Prints every action as it happens (set -x)
#   - Stops immediately on any error (set -e)
#
# USAGE:
#   1. cd into your repo root first: cd ~/Desktop/scrnaseq-latent-recovery
#   2. Copy this script into the repo root
#   3. Review it once (nothing here deletes any file's content, only moves)
#   4. Run: bash migrate_repo_structure.sh
#   5. Review `git status` and `git diff --stat --cached` before committing
#   6. Commit and push (final commands are printed at the end, not run
#      automatically -- you should look at the result first)
# =============================================================================

set -euo pipefail
set -x

# ---- Safety check: must be in a clean repo ----
if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: working tree is not clean. Commit or stash changes first."
  exit 1
fi

if [ ! -f "PROJECT_HANDOVER.md" ]; then
  echo "ERROR: this doesn't look like the repo root (PROJECT_HANDOVER.md not found)."
  echo "cd into the repo root first."
  exit 1
fi

# ---- Backup branch, just in case ----
BACKUP_BRANCH="pre_reorg_backup_$(date +%Y%m%d_%H%M%S)"
git branch "$BACKUP_BRANCH"
echo "Backup branch created: $BACKUP_BRANCH (you can always 'git checkout $BACKUP_BRANCH' to undo)"

# =============================================================================
# 1. Create the new top-level structure
# =============================================================================
mkdir -p code/01_simulation
mkdir -p code/02_real_data
mkdir -p code/03_preprocessing
mkdir -p code/04_metrics
mkdir -p code/05_phase_space
mkdir -p code/06_failure_modes
mkdir -p code/eda_checkpoint2
mkdir -p code/eda_checkpoint3
mkdir -p code/dev_history

mkdir -p results/step1_eda/tables
mkdir -p results/step1_eda/figures
mkdir -p results/step2_eda
mkdir -p results/step3_eda
mkdir -p results/step4_metrics
mkdir -p results/step5_phase_space/tables
mkdir -p results/step5_phase_space/figures
mkdir -p results/step6_failure_modes

# =============================================================================
# 2. Move R/ (Steps 1-2, EDA checkpoint 2) into code/
# =============================================================================
git mv R/01_simulation/*.R code/01_simulation/ 2>/dev/null || true
git mv R/02_real_data/*.R code/02_real_data/ 2>/dev/null || true
git mv R/02_real_data/*.sh code/02_real_data/ 2>/dev/null || true
git mv R/03_eda_checkpoint2/*.R code/eda_checkpoint2/ 2>/dev/null || true

# Remove now-empty R/ directory
rmdir R/01_simulation R/02_real_data R/03_eda_checkpoint2 R 2>/dev/null || true

# =============================================================================
# 3. Move scripts/ into code/, splitting final pipeline vs. dev_history
# =============================================================================

# --- eda_checkpoint3 ---
git mv scripts/eda_checkpoint3.R code/eda_checkpoint3/ 2>/dev/null || true

# --- Step 3 preprocessing (final versions only) ---
for f in step3_1_raw_pca_full.R step3_1_raw_pca_pilot.R step3_1_raw_recovery.R \
         step3_1_real_fix.R step3_2_libnorm_pca_full.R step3_3_log_pca_full.R \
         step3_4_shiftedlog_pca_full.R step3_5_sctransform_with_loadings.R \
         step3_6_glmpca_chunked_projection.R step3_6_glmpca_real.R \
         step3_6_glmpca_simulated.R step3_loadings_backfill.R \
         step3_10_manifest.R glmpca_projection_downstream_test.R \
         glmpca_projection_test.R; do
  [ -f "scripts/$f" ] && git mv "scripts/$f" code/03_preprocessing/
done
[ -d "scripts/gpu_glmpca" ] && git mv scripts/gpu_glmpca code/03_preprocessing/gpu_glmpca

# Superseded/rescue versions of Step 3.5 -> dev_history (kept as provenance)
for f in step3_5_sctransform_pca_full.R step3_5_sctransform_rescue.R \
         step3_5_sctransform_rescue2.R; do
  [ -f "scripts/$f" ] && git mv "scripts/$f" code/dev_history/
done

# --- Step 4 metrics ---
for f in step4_1to3_subspace_metrics.R step4_4to6_secondary_metrics.R \
         step4_9_compile_results.R; do
  [ -f "scripts/$f" ] && git mv "scripts/$f" code/04_metrics/
done

# --- Step 5 phase space ---
for f in step5_1_organize_by_axis.R step5_2to7_axis_plots.R \
         step5_8_critical_thresholds.R step5_9_phase_maps.R \
         step5_fix_scdesign3_dropout_labels.R; do
  [ -f "scripts/$f" ] && git mv "scripts/$f" code/05_phase_space/
done

# --- Step 6 failure modes (final pipeline versions) ---
for f in step6_1_technical_separation.R \
         step6_2_cluster_collapse.R step6_2_convert_null_control_replicates.R \
         step6_3_discover_null_control_files.R step6_3_phantom_clustering_v2.R \
         step6_4_variance_hijacking.R \
         step6_5a_pancreas_dpt.R step6_5b_simulated_log2fc_discovery.R \
         step6_5c_simulated_variance_ratio.R \
         step6_6_prebuild_real_reference_spaces.R step6_6_neighborhood_collapse.R \
         step6_7_subspace_rotation_slippage.R step6_7b_slippage_gradient_analysis.R \
         step6_8_cross_reference_boundaries.R; do
  [ -f "scripts/$f" ] && git mv "scripts/$f" code/06_failure_modes/
done
[ -f "run_step6_rerun.sh" ] && git mv run_step6_rerun.sh code/06_failure_modes/

# Superseded v1 of phantom clustering -> dev_history
[ -f "scripts/step6_3_phantom_clustering.R" ] && \
  git mv scripts/step6_3_phantom_clustering.R code/dev_history/

# --- Everything else remaining in scripts/ -> dev_history ---
# (kmeans_repro_test*, pilot_check*, diagnostic_schema_check*, debug_na,
#  alignment_check, fix_pbmc68k_true_group, and any script not explicitly
#  moved above)
if [ -d "scripts" ]; then
  for f in scripts/*.R scripts/*.py; do
    [ -f "$f" ] && git mv "$f" code/dev_history/ 2>/dev/null || true
  done
  rmdir scripts 2>/dev/null || true
fi

# =============================================================================
# 4. Reorganize results/ and figures/
# =============================================================================

# --- Step 1 EDA (results/tables + results/figures -> results/step1_eda) ---
[ -d "results/tables" ] && git mv results/tables/* results/step1_eda/tables/ 2>/dev/null || true
[ -d "results/figures" ] && git mv results/figures/* results/step1_eda/figures/ 2>/dev/null || true
rmdir results/tables results/figures 2>/dev/null || true

# --- Step 2 EDA (results/eda_checkpoint2 -> results/step2_eda) ---
[ -d "results/eda_checkpoint2" ] && git mv results/eda_checkpoint2/* results/step2_eda/ 2>/dev/null || true
rmdir results/eda_checkpoint2 2>/dev/null || true

# --- Loose files at results/ root -> step3 (GLM-PCA projection test outputs) ---
for f in glmpca_projection_downstream_result.rds glmpca_projection_test_result.rds; do
  [ -f "results/$f" ] && git mv "results/$f" results/step4_metrics/
done

# --- Step 3 EDA (top-level figures/eda_checkpoint_3 -> results/step3_eda) ---
[ -d "figures/eda_checkpoint_3" ] && git mv figures/eda_checkpoint_3/* results/step3_eda/ 2>/dev/null || true
rmdir figures/eda_checkpoint_3 2>/dev/null || true

# --- Step 5 figures (top-level figures/step5_* -> results/step5_phase_space/figures) ---
[ -d "figures/step5_9_phase_maps" ] && git mv figures/step5_9_phase_maps/* results/step5_phase_space/figures/ 2>/dev/null || true
[ -d "figures/step5_axis_sweeps" ] && git mv figures/step5_axis_sweeps/* results/step5_phase_space/figures/ 2>/dev/null || true
rmdir figures/step5_9_phase_maps figures/step5_axis_sweeps figures 2>/dev/null || true

# --- Step 6 results (data/processed/step6_*.csv -> results/step6_failure_modes) ---
# NOTE: these files currently live outside git (data/processed/ is gitignored).
# Copy them IN explicitly here since we now want them tracked.
for f in data/processed/step6_1_technical_separation.csv \
         data/processed/step6_2_cluster_collapse.csv \
         data/processed/step6_2_collapsed_pairs_detail.csv \
         data/processed/step6_3_phantom_clustering.csv \
         data/processed/step6_3_null_control_manifest.csv \
         data/processed/step6_4_variance_hijacking.csv \
         data/processed/step6_5a_pancreas_dpt.csv \
         data/processed/step6_5b_simulated_log2fc.csv \
         data/processed/step6_5c_simulated_variance_ratio.csv \
         data/processed/step6_6_neighborhood_collapse.csv \
         data/processed/step6_7_subspace_rotation_slippage.csv \
         data/processed/step6_7b_slippage_gradient_analysis.csv \
         data/processed/step6_8_cross_reference.csv \
         data/processed/step6_8_excluded_unreliable_thresholds.csv \
         data/processed/embedding_manifest.csv \
         data/processed/step4_master_results_table.csv \
         data/processed/step4_1to3_subspace_metrics.csv \
         data/processed/step4_4to6_secondary_metrics.csv \
         data/processed/step5_8_critical_thresholds.csv \
         data/processed/step5_axis_sweep_table.csv; do
  if [ -f "$f" ]; then
    dest="results/step6_failure_modes/$(basename "$f")"
    case "$(basename "$f")" in
      step4_*|embedding_manifest.csv) dest="results/step4_metrics/$(basename "$f")" ;;
      step5_*) dest="results/step5_phase_space/tables/$(basename "$f")" ;;
    esac
    cp "$f" "$dest"
    git add "$dest"
    echo "Copied and staged: $f -> $dest"
  else
    echo "WARNING: expected file not found, skipping: $f"
  fi
done

# =============================================================================
# 5. Add the new top-level files (must already exist in repo root -- see
#    accompanying instructions for LICENSE, RUNNING.md, and the EDA1 note)
# =============================================================================
[ -f "LICENSE" ] && git add LICENSE
[ -f "RUNNING.md" ] && git add RUNNING.md
[ -f "docs/eda_checkpoint1_reconstruction_note.md" ] && git add docs/eda_checkpoint1_reconstruction_note.md

# =============================================================================
# Done -- stop here and review before committing
# =============================================================================
set +x
echo ""
echo "============================================================"
echo "Migration complete. NOTHING HAS BEEN COMMITTED YET."
echo "============================================================"
echo "Next steps:"
echo "  1. Run: git status"
echo "  2. Run: git diff --cached --stat   (review what's staged)"
echo "  3. Manually check for any files git couldn't figure out"
echo "     (look for anything still sitting in old locations)"
echo "  4. When satisfied:"
echo "       git add -A"
echo "       git commit -m 'Reorganize repository structure for publication'"
echo "       git push origin main"
echo "  5. If anything looks wrong: git checkout $BACKUP_BRANCH"
echo "============================================================"
