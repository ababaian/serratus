#!/bin/bash
# worker.sh
#
# Base: 
# AMI :
#
# TODO: Add test to ensure docker is installed and can be run
#       on instance. Else Error exit.
# TODO: Ensure EC2 instance has correct IAM permissions prior
#       to committing to serratus-dl container pipeline

set -eux

# USAGE ===================================================
function usage {
  echo ""
  echo "Usage: worker.sh -c [serratus-dl:latest] -u [localhost:8000] -w 2 [OPTIONS]"
  echo ""
  echo "    Container Options"
  echo "    -c    container:version to run [serratus-dl:latest]"
  echo "    -A    optional string of arguments to pass to container"
  echo "          i.e  '-k s3://alternative_storage_bucket/' "
  echo ""
  echo "    Scheduler Information"
  echo "    -u    URL/IP to query for Scheduler webserver"
  echo ""
  echo "    Worker Options"
  echo "    -w    Number of simultanious jobs to process"
  echo ""
  echo "    -h     This usage message."
  exit 1
 }

# Default Parameters
CONTAINER="serratus-dl:latest"
CONTAINER_ARGS=''
SCHED=localhost:8000
WORKERS=2

# Read input
while getopts c:A:u:wh FLAG; do
  case $FLAG in
    # Input Options ---------
    c)
      CONTAINER=$(readlink -f $OPTARG)
      ;;
    A)
      CONTAINER_ARGS=$(readlink -f $OPTARG)
      ;;
    u)
      SCHED=$(readlink -f $OPTARG)
      ;;
    w)
      WORKERS=$(readlink -f $OPTARG)
      ;;
    h)  #show help ----------
      usage
      ;;
    \?) #unrecognized option - show help
      echo "Input parameter not recognized"
      usage
      ;;
  esac
done
shift $((OPTIND-1))


# FUNCTIONS ===============================================
# Worker process.  Run in a loop, grabbing data from the scheduler and
# processing it.
function wait_for_scheduler {
    while true; do
        curl -s "$SCHED/status" && break
        sleep 5
    done
}

function terminate_handler {
    echo "    $ACC_ID was terminated without completing. Reset status."
    echo "    In trap $(date -In)"
    # Tell server to reset this job to a "new" state, since we couldn't
    # finish processing it.
    curl -s "$SCHED/finish_split_job?job_id=$ACC_ID&status=new"
}

function main_loop {
    WORKER_ID=$1
    
    trap terminate_handler SIGUSR1
    # Note: The "& wait"s are important.  Without them, bash will wait for
    # the command to finish before executing its traps.  When we use "& wait",
    # the command will recieve the same trap (killing it), and then run our
    # trap handler, which tells the server our job failed.

    while true; do
        # Get a job
        echo "$WORKERID - Requesting job from Scheduler..."
        JOB_JSON=$(curl -s "$SCHED/start_split_job")
        ACTION=$(jq -r .action <(echo $JOB_JSON))

        case "$ACTION" in
            process)
                echo "  $WORKERID - Process State received."
                ;;
            wait)
                echo "  $WORKERID - Wait State received."
                sleep 10
                ;;
            shutdown)
                echo "  $WORKERID - Shutdown State received."
                exit 0
                ;;
            *)
                echo "  $WORKERID - ERROR: Unknown State received."
                exit 1
        esac

        # Parse SRA Accession ID
        ACC_ID=$(jq -r .acc_id <(echo $JOB_JSON))

        # Fire up the serratus-dl container
        # ENTRYPOINT: serratus-dl.sh
        # calls: scripts/serratus-dl/run_download.sh
        #        scripts/serratus-dl/run_split.sh
        #        scripts/serratus-dl/run_upload.sh

        echo "  $WORKERID - Run cmd: sudo docker run $CONTAINER -a -s $ACC_ID"
        echo ""
        sudo docker run $CONTAINER -a -s $ACC_ID & wait

        # Upload chunks
        # TODO: Figure out a way of counting uploaded fq-blocks from inside
        # the container and returnign those values here
        # N_fqblock_unpaired
        # N_fqblock_paired
        N=42
        echo "Uploading $N chunks. Yeehaw!"
        sleep 10 & wait

        # Tell the scheduler we're done
        echo "  $WORKERID - Job $ACC_ID is complete. Update scheduler."
        RESPONSE=$(curl -s "$SCHED/finish_split_job?job_id=$ACC_ID&status=split_done&N=$N")
        unset ACC_ID
    done
}

# SCRIPT ==================================================

# Boot ----------------------
wait_for_scheduler

# Work Functions ------------
for i in $(seq 1 "$WORKERS"); do
    main_loop "$i" & worker[i]=$!
done

sleep 1

# Spot Handling -------------
# Monitor AWS Cloudwatch for spot-termination signal
# if Spot termination signal detected, proceed with
# shutdown via SIGURS1 signal
SPOT_OK=true

# For a minor optimization; query the time of the
# spot-termination signal and shutdown in the last 10
# seconds to maximize chance the job finishes.

while $SPOT_OK; do
    sleep 5

    # Query for term-sig
    SPOT_TERM=$(jq -r .action <(echo $(curl http://169.254.169.254/latest/meta-data/spot/instance-action))

    if [ "$SPOT_TERM" -eq "stop" ]
    then
        echo "SPOT TERMINATION SIGNAL RECEIEVED."
        echo "Initiating shutdown procedures for all workers"

        for i in $(seq 1 "$WORKERS"); do
            kill -USR1 ${worker[i]} 2>/dev/null || true
        done
    fi
done


# End
echo "Kill jobs" $(date -In)
wait
echo "Wait done" $(date -In)

