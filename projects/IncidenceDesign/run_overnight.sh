#!/bin/bash
# =============================================================================
# run_overnight.sh
# Launch the tau-sweep MLE simulation and monitor progress every 10 minutes.
#
# Progress is tracked via checkpoint files (results/checkpoints/).
# Each completed (tau x incidence config) pair writes one checkpoint, so
# 0–25 checkpoints = 0–100% of total work (5 tau x 5 configs).
#
# What this script does:
#   1. Starts R in the background (nohup + caffeinate so it survives
#      screen lock, display sleep, and terminal close)
#   2. Saves the PID and log path to results/run_status.txt for reconnecting
#   3. Monitors progress in the foreground every 10 minutes
#      - Ctrl+C stops monitoring; R keeps running in the background
#      - Re-run ./run_overnight.sh next time to see current status
#        (it detects the existing process and skips re-launching)
#
# Usage:
#   ./run_overnight.sh            # first launch
#   ./run_overnight.sh            # re-run later to reconnect monitoring
#
# Manual progress check (no script needed):
#   ls results/checkpoints/ | wc -l   # 0-25 checkpoints = % done
#   tail -20 results/tau_sweep_run_*.log
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="${SCRIPT_DIR}/code"
RESULTS_DIR="${SCRIPT_DIR}/results"
CHECKPOINT_DIR="${RESULTS_DIR}/checkpoints"
STATUS_FILE="${RESULTS_DIR}/run_status.txt"

TOTAL_CHECKPOINTS=25   # 5 tau values x 5 incidence configs
INTERVAL_MINS=10
INTERVAL_SECS=$((INTERVAL_MINS * 60))

# ---------------------------------------------------------------------------
# Progress display
# ---------------------------------------------------------------------------
print_progress() {
  local start_epoch="$1"
  local now_epoch
  now_epoch=$(date +%s)
  local elapsed_secs=$(( now_epoch - start_epoch ))
  local elapsed_hrs=$(( elapsed_secs / 3600 ))
  local elapsed_mins=$(( (elapsed_secs % 3600) / 60 ))

  # Count completed checkpoints (= completed tau x config pairs)
  local n_done
  n_done=$(ls "${CHECKPOINT_DIR}"/*.rds 2>/dev/null | wc -l | tr -d ' ')
  n_done=${n_done:-0}

  local pct=$(( n_done * 100 / TOTAL_CHECKPOINTS ))

  # Progress bar (20 chars wide)
  local filled=$(( n_done * 20 / TOTAL_CHECKPOINTS ))
  local bar=""
  for (( i=0; i<filled; i++ )); do bar+="█"; done
  for (( i=filled; i<20; i++ )); do bar+="░"; done

  # Current tau from the log (last tau header seen)
  local log_file
  log_file=$(awk '{print $2}' "${STATUS_FILE}" 2>/dev/null | grep "LOG" | awk '{print $2}')
  local current_tau="—"
  if [ -n "${log_file}" ] && [ -f "${log_file}" ]; then
    local tau_line
    tau_line=$(grep "########## tau = " "${log_file}" 2>/dev/null | tail -1)
    if [ -n "${tau_line}" ]; then
      current_tau=$(echo "${tau_line}" | grep -oE '[0-9]+\.[0-9]+' | head -1)
    fi
  fi

  # Taus fully completed (all 5 configs done for that tau)
  local taus_done="?"
  if [ -d "${CHECKPOINT_DIR}" ]; then
    # Count distinct tau values that have 5 checkpoint files each
    taus_done=$(ls "${CHECKPOINT_DIR}"/*.rds 2>/dev/null \
      | sed 's/.*_tau//' | sed 's/\.rds//' \
      | sort | uniq -c \
      | awk '$1 == 5' | wc -l | tr -d ' ')
  fi

  # ETA estimate
  local eta_str="calculating..."
  if [ "${n_done}" -gt 0 ]; then
    local secs_per_unit=$(( elapsed_secs / n_done ))
    local remaining=$(( (TOTAL_CHECKPOINTS - n_done) * secs_per_unit ))
    local eta_hrs=$(( remaining / 3600 ))
    local eta_mins=$(( (remaining % 3600) / 60 ))
    eta_str="${eta_hrs}h ${eta_mins}m"
  fi

  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  printf "  [%s]  Elapsed: %dh %02dm\n" "$(date '+%H:%M')" "${elapsed_hrs}" "${elapsed_mins}"
  printf "  Progress:  [%s] %d/%d configs (%d%%)\n" "${bar}" "${n_done}" "${TOTAL_CHECKPOINTS}" "${pct}"
  printf "  Tau done:  %s/5 complete  |  Current tau: %s\n" "${taus_done}" "${current_tau}"
  printf "  ETA:       %s\n" "${eta_str}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}

# ---------------------------------------------------------------------------
# Detect existing run (reconnect mode)
# ---------------------------------------------------------------------------
if [ -f "${STATUS_FILE}" ]; then
  SAVED_PID=$(grep "^PID" "${STATUS_FILE}" | awk '{print $2}')
  SAVED_LOG=$(grep "^LOG" "${STATUS_FILE}" | awk '{print $2}')
  SAVED_START=$(grep "^START_EPOCH" "${STATUS_FILE}" | awk '{print $2}')

  if [ -n "${SAVED_PID}" ] && kill -0 "${SAVED_PID}" 2>/dev/null; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Reconnecting to existing run (PID ${SAVED_PID})"
    echo "  Log: ${SAVED_LOG}"
    echo "  (Ctrl+C to stop monitoring — R keeps running)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    print_progress "${SAVED_START}"
    while kill -0 "${SAVED_PID}" 2>/dev/null; do
      sleep "${INTERVAL_SECS}"
      kill -0 "${SAVED_PID}" 2>/dev/null && print_progress "${SAVED_START}"
    done
    echo ""
    echo "  Run complete! Final status:"
    print_progress "${SAVED_START}"
    rm -f "${STATUS_FILE}"
    exit 0
  else
    echo "Previous run (PID ${SAVED_PID}) is no longer active. Starting fresh."
    rm -f "${STATUS_FILE}"
  fi
fi

# ---------------------------------------------------------------------------
# Fresh launch
# ---------------------------------------------------------------------------
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${RESULTS_DIR}/tau_sweep_run_${TIMESTAMP}.log"
mkdir -p "${RESULTS_DIR}" "${CHECKPOINT_DIR}"

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Tau-Sweep Overnight Run"
echo "  Start: $(date)"
echo "  Log:   ${LOG_FILE}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Launching R (nohup + caffeinate -i -s)..."
echo "  Updates every ${INTERVAL_MINS} minutes. Ctrl+C stops monitoring;"
echo "  R continues running in background."
echo ""

nohup caffeinate -i -s \
  Rscript "${CODE_DIR}/05_run_simulation.R" \
  > "${LOG_FILE}" 2>&1 &

R_PID=$!
START_EPOCH=$(date +%s)

# Save state for reconnection
cat > "${STATUS_FILE}" <<EOF
PID ${R_PID}
LOG ${LOG_FILE}
START_EPOCH ${START_EPOCH}
EOF

echo "  PID: ${R_PID}"
echo "  To reconnect after closing this terminal: ./run_overnight.sh"
echo "  Quick check (no script): ls results/checkpoints/ | wc -l  (0-25)"
echo ""

# Brief pause so R can write its header before the first progress print
sleep 5

# Monitoring loop
print_progress "${START_EPOCH}"
while kill -0 "${R_PID}" 2>/dev/null; do
  sleep "${INTERVAL_SECS}"
  kill -0 "${R_PID}" 2>/dev/null && print_progress "${START_EPOCH}"
done

echo ""
echo "  ✓ Run complete! Final status:"
print_progress "${START_EPOCH}"
rm -f "${STATUS_FILE}"
