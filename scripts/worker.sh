#!/bin/bash
set -eux
# Worker process.  Run in a loop, grabbing data from the scheduler and
# processing it.
SCHED=localhost:8000
WORKERS=2

function wait_for_scheduler {
    while true; do
        curl -s "$SCHED/status" && break
        sleep 5
    done
}

function terminate_handler {
    echo "In trap $(date -In)"
    # Tell server to reset this job to a "new" state, since we couldn't
    # finish processing it.
    curl -s "$SCHED/finish_split_job?job_id=$ACC_ID&status=new"
}

function main_loop {
    trap terminate_handler SIGUSR1
    while true; do
        # Get a job
        JOB_JSON=$(curl -s "$SCHED/start_split_job")
        ACTION=$(jq -r .action <(echo $JOB_JSON))

        case "$ACTION" in
            process)
                ;;
            wait)
                sleep 10
                ;;
            shutdown)
                exit 0
                ;;
            *)
                echo ERROR Invalid State
                exit 1
        esac
        ACC_ID=$(jq -r .acc_id <(echo $JOB_JSON))

        # Download
        echo "Downloading doo dee dah!"
        sleep 1 & wait

        # Run split, but using flock to ensure only one thread at a time
        SPLIT_CMD=$(jq -r .split_cmd <(echo $JOB_JSON))
        echo "Running split command: $SPLIT_CMD"
        flock run_lock -c 'sleep 10' & wait

        # Upload chunks
        N=42
        echo "Uploading $N chunks. Yeehaw!"
        sleep 10 & wait

        # Tell the scheduler we're done
        RESPONSE=$(curl -s "$SCHED/finish_split_job?job_id=$ACC_ID&status=split_done&N=$N")
        unset ACC_ID
    done
}

wait_for_scheduler

# Start up a bunch of worker threads
for i in $(seq 1 "$WORKERS"); do
    main_loop & worker[i]=$!
done

# Testing testing
sleep 1

# while true; if $(curl metadata) == terminated; then
    for i in $(seq 1 "$WORKERS"); do
        kill -USR1 ${worker[i]} 2>/dev/null || true
    done
# done; fi; esac; END FOR; END TRANSACTION; </html>;

echo "Kill jobs" $(date -In)
wait
echo "Wait done" $(date -In)

