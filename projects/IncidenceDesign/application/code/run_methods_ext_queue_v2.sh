#!/bin/bash
# ============================================================
# Script: run_methods_ext_queue_v2.sh
# Purpose: Unattended queue for the 2026-10-21 reruns (author decisions of
#   2026-10-21): v2 checks, the A3 rerun with >= 2-college saturation regions and
#   t + df-corrected test counts (a3v2), its bit-for-bit reproduction check against
#   the 2026-10-20 A3 run, the long-chain Bayesian CAR rerun at rho = 0.5, the v2
#   summaries and V9. Unlike the 2026-10-20 queue, any failed step STOPS the queue
#   (plan: docs/plans/overnight_plan_2026-10-21.md).
# Author: Andrew Walther
# Created: 2026-10-21
# Dependencies: Rscript; methods_ext_*.R; application/tests/*
# ============================================================
# Run from the project root (keeps the machine awake while it runs):
#   caffeinate -i bash application/code/run_methods_ext_queue_v2.sh
# Status: application/results/methods_ext_2026-10/queue_v2.log ends with QUEUEDONE
# or QUEUEFAILED <step>. Each step's output: <results>/<step>.log.
# Compute: about 4 h (a3v2, 10 workers) + 4 h (a6bayeslong) on this machine.
set -u
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
R=application/results/methods_ext_2026-10
Q=$R/queue_v2.log
log() { echo "[$(date '+%F %T')] $*" >> "$Q"; }
stop_queue() { log "QUEUEFAILED $1 ($2)"; exit 1; }

# run <step> <success-pattern or -> <command...>: fails on non-zero exit, or when the
# success pattern (a line the step prints only when every task succeeded) is missing.
run() {
  name=$1; pat=$2; shift 2
  log "start $name"
  "$@" > "$R/${name}.log" 2>&1; rc=$?
  log "end $name (exit $rc)"
  [ $rc -eq 0 ] || stop_queue "$name" "exit $rc"
  if [ "$pat" != "-" ]; then grep -Eq "$pat" "$R/${name}.log" || stop_queue "$name" "missing '$pat'"; fi
}

log "queue v2 start"
run v2_tests "All v2 checks passed" Rscript application/tests/test_methods_ext_v2.R
run a3v2 "done: [0-9]+ OK, 0 failed" Rscript application/code/methods_ext_participation.R a3v2 10
run a3v2_reproduction "^PASS" Rscript application/tests/check_a3v2_reproduction.R a3 a3v2
run a6bayeslong "done: [0-9]+ fit rows, 0 failed tasks" Rscript application/code/methods_ext_sensitivity.R a6bayeslong 10
run summaries_v2 - Rscript application/code/methods_ext_summaries.R 10 a3v2 summaries_v2
run bayes_long_summary "Same data confirmed" Rscript application/code/methods_ext_bayes_long_summary.R
run v9_v2 "PASS" Rscript application/tests/verify_methods_ext_manifests.R
log "QUEUEDONE"
