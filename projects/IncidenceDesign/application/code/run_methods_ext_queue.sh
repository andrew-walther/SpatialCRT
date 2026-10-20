#!/bin/bash
# ============================================================
# Script: run_methods_ext_queue.sh
# Purpose: Unattended queue for the methods extensions after the A3/A4 run:
#   A3/A4 summaries, then A7, A5 Stage 2, A6 CAR/reverse, A6 Bayesian CAR,
#   the sensitivity summaries and V9. Each step logs separately; a failed step
#   is logged and the queue continues with the independent steps.
# Author: Andrew Walther
# Created: 2026-10-20
# Dependencies: Rscript; methods_ext_*.R
# ============================================================
# Run from the project root (keeps the machine awake while it runs):
#   nohup caffeinate -i bash application/code/run_methods_ext_queue.sh &
set -u
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
R=application/results/methods_ext_2026-10
Q=$R/queue.log
log() { echo "[$(date '+%F %T')] $*" >> "$Q"; }

log "waiting for A4 (ALLDONE in participation/a4_run.log)"
until grep -q ALLDONE "$R/participation/a4_run.log" 2>/dev/null; do sleep 60; done
log "A3/A4 finished"

run() { name=$1; shift; log "start $name"; "$@" > "$R/${name}.log" 2>&1; rc=$?; log "end $name (exit $rc)"; }
run summaries_a3_a4 Rscript application/code/methods_ext_summaries.R 10
run a7 Rscript application/code/methods_ext_sensitivity.R a7 10
run a5 Rscript application/code/methods_ext_sensitivity.R a5 10
run a6car Rscript application/code/methods_ext_sensitivity.R a6car 10
run sensitivity_summary_partial Rscript application/code/methods_ext_sensitivity_summary.R
run a6bayes Rscript application/code/methods_ext_sensitivity.R a6bayes 10
run sensitivity_summary Rscript application/code/methods_ext_sensitivity_summary.R
run v9 Rscript application/tests/verify_methods_ext_manifests.R
log "QUEUEDONE"
