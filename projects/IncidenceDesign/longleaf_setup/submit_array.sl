#!/bin/bash

# =============================================================================
# LONGLEAF JOB ARRAY TEMPLATE
# For running R simulations with different seeds/parameters in parallel
#
# Usage:
#   sbatch submit_array.sl
#
# Edit the CONFIGURATION section below before submitting.
# =============================================================================

# -----------------------------------------------------------------------------
# SLURM RESOURCE REQUEST
# Tip: start conservative, check seff <jobID> after first run, then tune.
# -----------------------------------------------------------------------------

#SBATCH --job-name=my_simulation        # name shown in squeue

#SBATCH --array=1-100                   # runs job 100 times (TASK_ID = 1..100)
                                        # change to e.g. 1-500, or 1-50%10
                                        # the %10 limits to 10 running at once

#SBATCH --ntasks=1                      # 1 task per array job (standard for R)
#SBATCH --cpus-per-task=1               # CPUs per task
                                        # increase to 4-8 if using parallel/future

#SBATCH --mem=8g                        # memory per task
                                        # start with 8g; check seff after first run

#SBATCH --time=12:00:00                 # max wall time per task (HH:MM:SS)
                                        # job is killed if it exceeds this

# Email notifications — replace with your actual UNC email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=onyen@email.unc.edu

# Log files — %A = array job ID, %a = task index
# Each task gets its own log file in the logs/ directory
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

# -----------------------------------------------------------------------------
# CONFIGURATION — edit these for your project
# -----------------------------------------------------------------------------

# Path to your R script (relative to where you run sbatch, or use full path)
R_SCRIPT="R/simulation.R"

# The SLURM task ID becomes the seed/parameter index passed to your R script
# Your R script reads this via: task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
TASK_ID=${SLURM_ARRAY_TASK_ID}

# -----------------------------------------------------------------------------
# ENVIRONMENT SETUP
# -----------------------------------------------------------------------------

module purge
module add r/4.4.0

# Print job info to log (helpful for debugging)
echo "=============================="
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "=============================="

# -----------------------------------------------------------------------------
# RUN
# -----------------------------------------------------------------------------

Rscript "${R_SCRIPT}" "${TASK_ID}"

echo "Finished at: $(date)"
