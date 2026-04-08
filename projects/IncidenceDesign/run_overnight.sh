#!/bin/bash
# =============================================================================
# run_overnight.sh
# Launch the tau-sweep MLE simulation locally, running overnight.
#
# What this does:
#   - Changes to the code directory so 05_run_simulation.R can find its modules
#   - Redirects all R output (stdout + stderr) to a timestamped log file
#   - Uses nohup so the job survives terminal/SSH disconnect
#   - Uses caffeinate so macOS does not sleep and kill the process
#   - Prints the PID so you can monitor or kill it if needed
#
# Expected runtime: ~14-15 hours (5 cores x 5 incidence configs x 5 tau values)
#
# Usage:
#   chmod +x run_overnight.sh    # first time only
#   ./run_overnight.sh           # launch from ANY working directory
#
# Monitor progress:
#   tail -f results/tau_sweep_run_<timestamp>.log
#
# Check if still running:
#   ps aux | grep 05_run_simulation
#   # or use the PID printed at launch:
#   ps -p <PID>
#
# Kill the job (if needed):
#   kill <PID>
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="${SCRIPT_DIR}/code"
RESULTS_DIR="${SCRIPT_DIR}/results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${RESULTS_DIR}/tau_sweep_run_${TIMESTAMP}.log"

# Ensure results dir exists (it should already)
mkdir -p "${RESULTS_DIR}"

echo "=============================="
echo "  Tau-Sweep Overnight Run"
echo "=============================="
echo "Script dir:  ${SCRIPT_DIR}"
echo "Log file:    ${LOG_FILE}"
echo "Start time:  $(date)"
echo ""
echo "Launching R process..."

# Launch:
#   - nohup: survive terminal close
#   - caffeinate -i: prevent macOS idle sleep for duration of process (-w PID added below)
#   - Rscript: run the simulation script
#   - All output (stdout + stderr) goes to the log file
nohup caffeinate -i \
  Rscript "${CODE_DIR}/05_run_simulation.R" \
  > "${LOG_FILE}" 2>&1 &

R_PID=$!
echo "R process PID: ${R_PID}"
echo ""
echo "To monitor progress:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "To check if still running:"
echo "  ps -p ${R_PID}"
echo ""
echo "To stop the run:"
echo "  kill ${R_PID}"
echo ""
echo "Expected completion: ~14-15 hours from now (~$(date -v +15H '+%H:%M tomorrow'))"
echo "=============================="
