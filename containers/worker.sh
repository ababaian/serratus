#!/usr/bin/bash
set -eu

# Wrapper script to serratus-{dl,align-merge}.  This script provides looping,
# multi-threading and also checks for spot termination.  If nodes are
# terminated by Spot, then it shuts everything down and tells the scheduler.

if [ ! -x "$2" ]; then
    echo "usage: worker.sh <split|align|merge> <script.sh> [args]"
    exit 1
fi
TYPE="$1"; shift

#if [ "$TYPE" = "merge" ]; then
#    # protect from termination
#    touch running.merge
#fi

if [ -z "$SCHEDULER" ]; then
    echo Please set SCHEDULER environment variable.
    exit 1
fi

# Run with nproc workers by default.  We can probably improve CPU usage
# by running 2*nproc, at the expense of disk space...  which is a bit
# limiting at the moment.
WORKERS=${WORKERS:-$(nproc)}
retry_count=0

function terminate_handler {
    # Tell server to reset this job to a "new" state,
    # since we couldn't finish processing it.
    curl -s -X POST "$SCHEDULER/jobs/$TYPE/$JOB_ID&state=new"
    
    echo "    $JOB_ID was terminated without completing. Reset status."
    echo "    In trap $(date -In)"

}

function main_loop {
    trap terminate_handler SIGUSR1
    # Note: The "& wait"s are important.  Without them, bash will wait for
    # the command to finish before executing its traps.  When we use "& wait",
    # the command will recieve the same trap (killing it), and then run our
    # trap handler, which tells the server our job failed.
    WORKER_ID="$INSTANCE_ID-$1"
    NWORKER="$1"
    shift

    # To off-set AWS API calls, offset all workers coming online
    # at once. Sleep for 0.0-9.9 seconds once.
    if [ "$FIRST" = "true" ]
    then
        sleep $[ ( $RANDOM % 9 ) ].$[ ( $RANDOM % 9 ) ]
        FIRST='false'
    fi


    while true; do
        echo "$WORKER_ID - Requesting job from Scheduler..."

        if [ -f "$BASEDIR/.shutdown-lock" ]
        then
            # Shutdown on another worker/thread is initiated.
            # Do not request work. Shutdown imminent
            sleep 300s
        fi

        JOB_JSON=$(curl -fs -X POST "$SCHEDULER/jobs/$TYPE/?worker_id=$WORKER_ID" || true)

        if [ "$TYPE" = align ]; then
            JOB_ID=$(echo $JOB_JSON | jq -r .block_id)
        else
            JOB_ID=$(echo $JOB_JSON | jq -r .acc_id)
        fi

        if [ -n "$JOB_JSON" ]; then
            ACTION=$(echo $JOB_JSON | jq -r .action)
        else
            ACTION=retry
        fi

        # Maximum failed retries before
        # self-termination is initiated
        if [ $retry_count -gt 5 ]
        then
            ACTION=shutdown
        fi

        case "$ACTION" in
          process)
            echo "  $WORKER_ID - Process State received.  Running $@"
            retry_count=0 # reset retry counter

            export JOB_JSON WORKER_ID

            # Worker punch-in
            # offset by worker # seconds to not collide dl queries
            if [ ! -f "running.$NWORKER" ]; then
                touch running.$NWORKER

                if [ ! -f $BASEDIR/scale.in.pro ]; then
                  echo "  Enable SCALE-IN protection"
                  touch $BASEDIR/scale.in.pro

                  # Turn ON scale-in protection
                  aws autoscaling set-instance-protection \
                    --region us-east-1 \
                    --instance-ids $INSTANCE_ID \
                    --auto-scaling-group-name $ASG_NAME \
                    --protected-from-scale-in & wait
                    
                fi
            fi

            # Run the target script.
            "$@" & wait

            # Unset job ID to prevent incorrect terminations
            unset JOB_ID
            ;;
          wait)
            echo "  $WORKER_ID - Wait State received."

            # Worker punch-out
            if [ -f "running.$NWORKER" ]; then
                rm -f running.$NWORKER
            fi

            # When all worker's are punched-out
            # remove scale-in protection
            if ls running* 1> /dev/null 2>&1
            then
                ## Other worker is not checked out
                retry_count=0
            elif [ -f "$BASEDIR/scale.in.pro" ]
            then
                echo "  Removing SCALE-IN protection"
                # Turn off scale-in protection
                aws autoscaling set-instance-protection \
                  --region us-east-1 \
                  --instance-ids $INSTANCE_ID \
                  --auto-scaling-group-name $ASG_NAME \
                  --no-protected-from-scale-in & wait

                rm -f $BASEDIR/scale.in.pro
            fi

            # Add to retry counter
            ((retry_count=retry_count+1))

            # wait cycle
            sleep 10 & wait

            continue
            ;;
          retry)
            echo "  $WORKER_ID - Scheduler not reached.  Waiting"

            if ls running* 1> /dev/null 2>&1; then
                ## Other worker is not checked out
                retry_count=0
            else
                # Add to retry counter
                ((retry_count=retry_count+1))
            fi

            sleep 10 & wait
            continue
            ;;
          shutdown)
            echo "  $WORKER_ID - Shutdown State received."

            # This is not threadsafe!  For now let's just put it in a critical section.
            # Using a recipe from "man flock" which appears to work.
            (
                flock 200

                echo "  Shutting down instance"
                # TODO: change to shutdown (see below)
                aws ec2 terminate-instances \
                 --region us-east-1 \
                 --instance-ids $INSTANCE_ID

                sleep 300
                false
                exit 0

            ) 200> "$BASEDIR/.shutdown-lock"
           
            ;;
          *)        echo "  $WORKER_ID - ERROR: Unknown State received."
            echo "  $WORKER_ID - ERROR: Unknown State received."
            # Don't exit; we just got one invalid request.
        esac
    done
}

echo "=========================================="
echo "              SERRATUS  INIT              "
echo "=========================================="
cd $BASEDIR

# AWS Internal check ------------
INSTANCE_ID=$(curl -s http://169.254.169.254/latest/meta-data/instance-id)

# Check AWS Credential
aws s3 cp s3://serratus-public/var/aws-test-token.jpg ./ || true

if [ ! -f ./aws-test-token.jpg ];
then
    # Internal credential error - kill
    # Paradox, can't terminate without credentials
    # --> solution: change default shutdown behavior
    #     to "terminate" in terraform for all ASG
    #

    sudo shutdown

    #aws ec2 terminate-instances \
    # --region us-east-1 \
    # --instance-ids $INSTANCE_ID

    false
    exit 0
fi

# Retrieve ASG name
ASG_NAME=$(aws ec2 describe-tags --filters "Name=resource-id,Values=$INSTANCE_ID" --region us-east-1 | jq -r '.Tags[] | select(.["Key"] | contains("aws:autoscaling:groupName")) | .Value')
FIRST='true'
export INSTANCE_ID ASG_NAME FIRST

# ----------------------------

# Fire up main loop (SRA downloader)
echo "Creating $WORKERS worker processes"
for i in $(seq 1 "$WORKERS"); do
    main_loop "$i" "$@" & worker[i]=$!
done

function kill_workers {
    for i in $(seq 1 "$WORKERS"); do
        kill -USR1 ${worker[i]} 2>/dev/null || true
    done
    
    exit 0
}

# Send signal if docker is shutting down.
trap kill_workers TERM

# Spot Operations ===============================
# Monitor AWS Cloudwatch for spot-termination signal
# if Spot termination signal detected, proceed with
# shutdown via SIGURS1 signal

METADATA=http://169.254.169.254/latest/meta-data
while true; do
    # Note: this URL returns an HTML 404 page when there is no action.  Use
    # "curl -f" to mitigate that.
    INSTANCE_ACTION=$(curl -fs $METADATA/spot/instance-action | jq -r .action)
    if [ "$INSTANCE_ACTION" == "terminate" ]; then
        echo "SPOT TERMINATION SIGNAL RECEIEVED."
        echo "Initiating shutdown procedures for all workers"

        kill_workers
    fi

    sleep 5
done
