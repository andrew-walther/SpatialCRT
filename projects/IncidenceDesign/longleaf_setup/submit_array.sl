#!/bin/bash

# =============================================================================
# submit_array.sl — SLURM job array for IncidenceDesign MLE tau-sweep simulation
#
# Runs 12,800 scenarios in parallel (one per array task).
# 12,800 = 5 tau values x 2,560 base scenarios (25 x 10 resamples each).
# Each task takes ~1-5 min with default resamples (25 x 10).
#
# Usage:
#   sbatch submit_array.sl                    # submit all 12,800 tasks
#   sbatch --array=1-10 submit_array.sl       # dry run with first 10 tasks only
#   sbatch --array=1-1 submit_array.sl        # test with task 1 only
#   sbatch --array=5,12,47 submit_array.sl    # resubmit specific failed tasks
# =============================================================================

# --- SLURM Resource Request ---
#SBATCH --job-name=spatialcrt_tau_sweep
#SBATCH --array=1-12800%100             # 12,800 tasks, max 100 concurrent
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g                        # tune with seff after first run
#SBATCH --time=00:30:00                 # Phase 2: increase to 01:00:00

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=awalther@email.unc.edu

# Log files (%A = array job ID, %a = task index)
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

# --- Environment ---
module purge
module add r/4.4.0

# --- Job Info ---
echo "=============================="
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "=============================="

# --- Run ---
Rscript simulation.R ${SLURM_ARRAY_TASK_ID}

echo "Finished at: $(date)"
