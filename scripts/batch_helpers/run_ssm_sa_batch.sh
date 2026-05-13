#!/usr/bin/env bash
# Run all 60 H_* SSM-SA instances sequentially using detached screen sessions.
# Each instance is launched only after the previous screen session exits.
#
# Parameters are set to irace Config 4 (best rank):
#   T0=4.5479  Tf=0.0373  Nit=150
#
# IMPORTANT (remote server): run this script itself inside a persistent session
# so that it survives SSH disconnects, e.g.:
#
#   screen -S ssm_sa_batch -dmL -Logfile logs_ssm_sa/ssm_sa_batch.log \
#       bash scripts/run_ssm_sa_batch.sh
#
# or:
#   nohup bash scripts/run_ssm_sa_batch.sh > logs_ssm_sa/ssm_sa_batch.log 2>&1 &
#
# Usage:
#   bash scripts/run_ssm_sa_batch.sh [instances_file]
#   bash scripts/run_ssm_sa_batch.sh scripts/list_ssm_sa_all.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

INSTANCES_FILE=""
BETA_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --beta=*) BETA_ARG="$1"; shift ;;
        --beta)   BETA_ARG="--beta=$2"; shift 2 ;;
        -*)       echo "Unknown option: $1" >&2; exit 1 ;;
        *)        INSTANCES_FILE="$1"; shift ;;
    esac
done
INSTANCES_FILE="${INSTANCES_FILE:-"$SCRIPT_DIR/list_ssm_sa_all.txt"}"}

LOG_DIR="$PROJECT_DIR/logs_ssm_sa"
POLL_INTERVAL=60   # seconds between screen-session checks

mkdir -p "$LOG_DIR"

echo "=== SSM-SA batch run ==="
echo "Instances file  : $INSTANCES_FILE"
echo "Log directory   : $LOG_DIR"
echo "Project dir     : $PROJECT_DIR"
echo "Poll interval   : ${POLL_INTERVAL}s"
echo "Parameters      : T0=4.5479  Tf=0.0373  Nit=150  (irace Config 4)"
[[ -n "$BETA_ARG" ]] && echo "Beta override   : $BETA_ARG"
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
        "$SCRIPT_DIR/run_ssm_sa_cli.jl" "--instance=$instance" ${BETA_ARG:+$BETA_ARG}

    # Poll until the screen session disappears
    while screen -list | grep -qF ".$session"; do
        sleep "$POLL_INTERVAL"
    done

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed : $session"
    echo ""
done < "$INSTANCES_FILE"

echo "=== All instances completed at $(date) ==="
