#!/usr/bin/env bash
# Run MILP instances sequentially using detached screen sessions.
# Each instance is launched only after the previous screen session exits.
#
# Expected runtime per instance: ~30000 s (~8.5 h). Polling interval: 300 s.
#
# IMPORTANT (remote server): run this script itself inside a persistent session
# so that it survives SSH disconnects, e.g.:
#
#   screen -S milp_batch -dmL -Logfile logs_b3/milp_batch.log \
#       bash scripts/run_milp_batch.sh
#
# or:
#   nohup bash scripts/run_milp_batch.sh > logs_b3/milp_batch.log 2>&1 &
#
# Usage:
#   bash scripts/run_milp_batch.sh [instances_file]
#   bash scripts/run_milp_batch.sh scripts/lowerbound_beta_3_list.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
INSTANCES_FILE="${1:-"$SCRIPT_DIR/list_complete_beta3.txt"}"
LOG_DIR="$PROJECT_DIR/logs_beta3"
POLL_INTERVAL=30   # seconds between screen-session checks (~8.5 h runs)

mkdir -p "$LOG_DIR"

echo "=== MILP batch run ==="
echo "Instances file  : $INSTANCES_FILE"
echo "Log directory   : $LOG_DIR"
echo "Project dir     : $PROJECT_DIR"
echo "Poll interval   : ${POLL_INTERVAL}s"
echo "Started at      : $(date)"
echo ""

while IFS= read -r instance || [[ -n "$instance" ]]; do
    instance="${instance%$'\r'}"

    # Skip blank lines and comment lines
    if [[ -z "$instance" || "$instance" == \#* ]]; then
        continue
    fi

    session="${instance%.jl}"
    logfile="$LOG_DIR/${session}.log"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting  : $session"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file  : $logfile"

    # screen -list shows sessions as "PID.SESSION_NAME"; match on ".SESSION_NAME"
    screen -S "$session" -dmL -Logfile "$logfile" \
        julia --project="$PROJECT_DIR" \
        "$SCRIPT_DIR/run_milp_cli.jl" "--instance=$instance"

    # Poll until the screen session disappears
    while screen -list | grep -qF ".$session"; do
        sleep "$POLL_INTERVAL"
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished  : $session"
    echo ""

done < "$INSTANCES_FILE"

echo "=== All instances completed at $(date) ==="
