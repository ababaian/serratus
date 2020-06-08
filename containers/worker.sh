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

function main_loop {
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
        if [ -f "$BASEDIR/.shutdown-lock" ]
        then
            # Shutdown on another worker/thread is initiated.
            # Do not request work. Shutdown imminent
            sleep 300s
        fi

        JOB_JSON=$(curl -fs -X POST "$SCHEDULER/jobs/$TYPE/?worker_id=$WORKER_ID" || true)

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
            fail_count=$retry_count  # temporary holder for retry_counter
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
                    --protected-from-scale-in
                    
                fi
            fi

            # If `serratus-align.sh` fails
            # it will `touch $RUN_FAIL`
            RUN_FAIL="$BASEDIR"/"$NWORKER".fail
            export RUN_FAIL

            # Run the target script.
            "$@"

            # If run-failed. Count the failure as a 
            # 'retry_count' to initiate shut-down
            #
            if [ -e "$RUN_FAIL" ]; then
                ((retry_count=fail_count+1))
                rm $RUN_FAIL
            fi
            ;;
          wait)
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
                  --no-protected-from-scale-in

                rm -f $BASEDIR/scale.in.pro
            fi

            # Add to retry counter
            ((retry_count=retry_count+1))

            # wait cycle
            sleep 10

            continue
            ;;
          retry)
            if ls running* 1> /dev/null 2>&1; then
                ## Other worker is not checked out
                retry_count=0
            else
                # Add to retry counter
                ((retry_count=retry_count+1))
            fi

            sleep 10
            continue
            ;;
          shutdown)
            (
                flock 200

                echo "  Shutting down instance"
                # TODO: change to shutdown (see below)
                aws ec2 terminate-instances \
                 --region us-east-1 \
                 --instance-ids $INSTANCE_ID

                sleep 300

                # TODO: Add a redundancy for shutdown
                #       to work form inside the container
                #
                # Secondary back-up -- shutdown instance
                # (set to "stopped" state" if terminate fails)
                # yum -y install sudo shadow-utils util-linux
                # sudo shutdown -h now
                # sleep 300

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

    sudo shutdown
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

# :)
