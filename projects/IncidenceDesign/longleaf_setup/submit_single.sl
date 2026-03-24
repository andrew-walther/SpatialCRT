#!/bin/bash

# =============================================================================
# LONGLEAF SINGLE JOB TEMPLATE
# For one long-running R simulation (not an array).
# Use this when you have a single script that runs start-to-finish.
# Use submit_array.sl instead when running many seeds/parameters.
# =============================================================================

#SBATCH --job-name=my_single_sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1               # increase if using parallel/future in R
#SBATCH --mem=16g
#SBATCH --time=2-00:00:00               # 2 days; format: D-HH:MM:SS

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=onyen@email.unc.edu

#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

# -----------------------------------------------------------------------------
module purge
module add r/4.4.0

echo "Job ID:   ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Start:    $(date)"

Rscript R/simulation.R

echo "End: $(date)"
