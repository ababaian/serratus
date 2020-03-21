#!/usr/bin/bash
set -eux
if [ ! -x "$2" ]; then
    echo "usage: worker.sh <split|align|merge> <script.sh> [args]"
    exit 1
fi

TYPE="$1"; shift

# Default base directory is /home/serratus
function wait_for_scheduler {
    while true; do
        if [ "$(curl -s "$SCHED/status" | jq -r .status)" == "up" ]; then
            break
        else
            sleep 5
        fi
    done
}

function terminate_handler {
    echo "    $ACC_ID was terminated without completing. Reset status."
    echo "    In trap $(date -In)"
    # Tell server to reset this job to a "new" state, since we couldn't
    # finish processing it.
    curl -s -X POST "$SCHED/jobs/$TYPE/$ACC_ID&state=terminated"
}

function main_loop {
    trap terminate_handler SIGUSR1
    # Note: The "& wait"s are important.  Without them, bash will wait for
    # the command to finish before executing its traps.  When we use "& wait",
    # the command will recieve the same trap (killing it), and then run our
    # trap handler, which tells the server our job failed.
    # what if we & wait on the whole main_loop function instead? Test this.
    WORKERID=$1

    # TODO: Wrap job query into self-contained function?
    while true; do
        echo "$WORKERID - Requesting job from Scheduler..."
        JOB_JSON=$(curl -s -X POST "$SCHED/jobs/$TYPE/")
        ACTION=$(echo $JOB_JSON | jq -r .action)

        case "$ACTION" in
          process)
            echo "  $WORKERID - Process State received.  Running $@"
            export JOB_JSON
            "$@" & wait
            ;;
          wait)
            echo "  $WORKERID - Wait State received."
            sleep 3
            continue
            ;;
          shutdown)
            echo "  $WORKERID - Shutdown State received."
            exit 0
            ;;
          *)        echo "  $WORKERID - ERROR: Unknown State received."
            exit 1
            echo "  $WORKERID - ERROR: Unknown State received."
            exit 1
        esac
    done
}

echo "=========================================="
echo "                SERRATUS                  "
echo "=========================================="
cd $BASEDIR
wait_for_scheduler

# Fire up main loop (SRA downloader)
for i in $(seq 1 "$WORKERS"); do
    main_loop "$i" & worker[i]=$!
done

# Spot Operations ===============================
# Monitor AWS Cloudwatch for spot-termination signal
# if Spot termination signal detected, proceed with
# shutdown via SIGURS1 signal

# TODO: For a minor optimization; query the time of the
# spot-termination signal and shutdown in the last 10
# seconds to maximize chance the job finishes.

METADATA=http://169.254.169.254/latest/meta-data
while true; do
    # Note: this URL returns an HTML 404 page when there is no action.  Use
    # "curl -f" to mitigate that.
    INSTANCE_ACTION=$(curl -fs $METADATA/spot/instance-action | jq -r .action)
    if [ "$INSTANCE_ACTION" == "terminate" ]; then
        echo "SPOT TERMINATION SIGNAL RECEIEVED."
        echo "Initiating shutdown procedures for all workers"

        for i in $(seq 1 "$WORKERS"); do
            kill -USR1 ${worker[i]} 2>/dev/null || true
        done
        break
    fi

    sleep 60 # TODO Change this back to 5.
done
