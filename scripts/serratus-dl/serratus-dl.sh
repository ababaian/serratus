#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-dl
# =====================================
set -ex
#
# Base: serratus-Downloader (>0.1.1)
# Amazon Linux 2 with Docker
# AMI : aami-0fdf24f2ce3c33243
# login: ec2-user@<ipv4>
# base: 9 Gb
# Container: serratus-dl:0.1
#

# Test cmd: ./serratus-dl.sh -s SRR11166696
# TODO: Consider switching to an external definition for RUNID
# such that downstream scripts can easily access the 
# $RUNID/ folder and it's files.
#
# TODO: run_split.sh script assumes a single fastq file (-f) is
# single-end reads and not interleaved or mixed paire-end reads.
# Adapt script to deal with these edge cases
#
PIPE_VERSION="0.1"
AMI_VERSION='ami-0fdf24f2ce3c33243'
CONTAINER_VERSION='serratus-dl:0.1'

# Usage
function usage {
  echo ""
  echo "Usage: docker exec serratus-dl -u [localhost:8000] -k <s3://bucket/fq-blocks> [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Fields"
  echo "    -u    <IP>:<PORT> to query for scheduler webserver [localhost:8000]"
  echo "    -k    S3 bucket path for fq-blocks [s3://serratus-public/fq-blocks]"
  echo ""
  echo "    Scheduler Information"
  echo "    -w    (N/A) Number of parallel download operations to run [1]"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo ""
  echo "    SRA Accession for sratools"
  echo "    -s    (N/A) SRA Accession for direct-to-pipeline access"
  echo ""
  echo "    AWS / S3 Bucket parameters"
  echo "    -a    Flag. Imports AWS IAM from this ec2-instance for container"
  echo "          EC2 instance must have been launched with correct IAM"
  echo "          (No alternative yet, hard set to TRUE)"
  echo ""
  echo "    Arguments as string to pass to downstream scripts"
  echo "    -D    String of arguments to pass to 'run_download.sh <ARG>'"
  echo "    -P    String of arguments to pass to 'run_split.sh <ARG>'"
  echo "    -U    String of arguments to pass to 'run_upload.sh <ARG>'"
  echo ""
  echo "    Output options"
  echo "    -d    Base directory in container [/home/serratus/]"
  echo "    -o    <output_filename_prefix> [Defaults to SRA_ACCESSION]"
  echo ""
  echo "    Outputs Uploaded to s3: "
  echo "          <output_prefix>.fq.xxxx ... <output_prefix>.fq.yyyy"
  echo ""
  echo "ex: docker exec serratus-dl -u localhost:8000"
  exit 1
}

# PARSE INPUT =============================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Scheduler / Container Parameters
SCHED='localhost:8000'
WORKERS='1'

# SRA Accession -S
SRA=''
THREADS='1'

# AWS Options -ak
AWS_CONFIG='TRUE'
S3_BUCKET='s3://serratus-public/fq-blocks'

# Script Arguments -DPU
DL_ARGS=''
SPLIT_ARGS=''
UL_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts u:w:n:s:ak:D:P:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    u)
      SCHED=$OPTARG
      ;;
    w)
      WORKERS=$OPTARG
      ;;
    n)
      THREADS=$OPTARG
      ;;
    s)
      SRA=$OPTARG
      ;;
    a)
      AWS_CONFIG='TRUE'
      ;;

    k)
      S3_BUCKET=$OPTARG
      ;;
    # SCRIPT ARGUMENTS -------
    D)
      DL_ARGS=$OPTARG
      ;;
    P)
      SPLIT_ARGS=$OPTARG
      ;;
    U)
      UL_ARGS='TRUE'
      ;;
    # output options -------
    d)
      BASEDIR=$OPTARG
      ;;
    o)
      OUTNAME=$OPTARG
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

# Check inputs --------------
if [[ ( -z "$OUTNAME" ) ]]; then OUTNAME="$SRA"; fi

# FUNCTIONS ===============================================
# Worker process.  Run in a loop, grabbing data from the scheduler and
# processing it.

function boot_procedure {
  # Pass credentials to sra-toolkit via vdb-config
  # Current version requires a manual
  # vdb-config -i initialization
  # ugly hack is to copy a blank config file and bypass this
  # TODO: Move vdb-config -i hack to Dockerfile/installation
  mkdir -p /root/.ncbi
  vdb-config --report-cloud-identity yes

  ## Download AWS S3 test token
  #aws s3 cp s3://serratus-public/aws-test-token.jpg $WORKDIR/aws-test-token.jpg

  #if [ ! -f "$WORKDIR/aws-test-token.jpg" ]
  #then
  #  echo "    ERROR: AWS Test did not download"
  #  echo "    Ensure the EC2 instance has correct IAM permissions"
  #  echo "      - requires S3 Read/Write"
  #  usage
  #fi

  #echo '  ...authentication token was download successfully'
  #echo ""

}

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
    curl -s "$SCHED/jobs/split/$ACC_ID&status=new"
}

function main_loop {
    trap terminate_handler SIGUSR1
    # Note: The "& wait"s are important.  Without them, bash will wait for
    # the command to finish before executing its traps.  When we use "& wait",
    # the command will recieve the same trap (killing it), and then run our
    # trap handler, which tells the server our job failed.
    # what if we & wait on the whole main_loop function instead? Test this.
    WORKER_ID=$1

    # TODO: Wrap job query into self-contained function?
    while true; do
        echo "$WORKERID - Requesting job from Scheduler..."
        JOB_JSON=$(curl -s -X POST "$SCHED/jobs/split/")
        ACTION=$(echo $JOB_JSON | jq -r .action)

        case "$ACTION" in
          process)
            echo "  $WORKERID - Process State received."
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

        # Parse SRA Accession ID
        ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)
        SRA_RUN=$(echo $JOB_JSON | jq -r .sra_run_id.Run)

        # TODO: Allow the scheduler/main data-table to have arugments
        # which will be passed on to the downloader scripts
        
        # TODO: If we're going to lock the system so the split commands
        # in each container don't collide on the CPU then we need to mount
        # a folder from filesystem into the container and `flock` a file
        # on the shared mount. For now use 1 worker per download EC2?

        # Run job --------------------------------------------------
        # ----------------------------------------------------------
        # Generate random alpha-numeric for run-id
        RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )
        WORKDIR=$BASEDIR/$RUNID
        mkdir -p $WORKDIR; cd $WORKDIR

        echo "============================"
        echo "    serratus-dl Pipeline    "
        echo "============================"
        echo " date:      $(date)"
        echo " version:   $PIPE_VERSION"
        echo " ami:       $AMI_VERSION"
        echo " container: $CONTAINER_VERSION"
        echo " worker-id: $WORKER_ID"
        echo " run-id:    $RUNID"
        echo " sra:       $SRA_RUN"
        echo " S3 url:    $S3_BUCKET"
        echo ""

        # RUN DOWNLOAD ==============================================
        echo "  Running -- run_download.sh --"
        echo "  $BASEDIR/scripts/run_download.sh -s $SRA_RUN $DL_ARGS"

        ./run_download.sh -s $SRA_RUN -p $THREADS $DL_ARGS & wait

        echo ''

        # Detect downloaded fastq files for split logic
        FQ0=$(ls *_0.fastq 2>/dev/null || true)
        FQ1=$(ls *_1.fastq 2>/dev/null || true)
        FQ2=$(ls *_2.fastq 2>/dev/null || true)

        if [[ ( -s $FQ1 && -n $FQ1 ) && ( -s $FQ2 && -n $FQ2 ) ]]
        then
          paired_exists=true
          echo "  Paired-end reads detected"
        else
          paired_exists=false
          echo "  Paired-end reads not-detected"
        fi

        if [[ ( -s $FQ0 && -n $FQ0 ) ]]
        then
          unpaired_exists=true
          echo "  Unpaired reads detected"
        else
          unpaired_exists=false
          echo "  unpaired reads not-detected"
        fi

        if [[ "$paired_exists" = true && "$unpaired_exists" = true ]]
        then
          # Both paired-end and unpaired exist
          # use only paired-end data
          echo "   WARNING: Paired and Unpaired data detected"
          echo "            Using only Paired-End Reads"
          unpaired_exists=FALSE
        fi

        # RUN SPLIT ===============================================
        # Add FQ0 vs. FQ1+FQ2 logic here
        echo "  Running -- run_split.sh --"

        if [[ "$paired_exists" = true ]]
        then
          echo "  .$BASEDIR/scripts/run_split.sh -o $OUTNAME -p $THREADS $SPLIT_ARGS"
          bash $BASEDIR/scripts/run_split.sh -1 $FQ1 -2 $FQ2 -o $SRA_RUN -p $THREADS $SPLIT_ARGS & wait

        elif [[ "$paired_exists" = true ]]
        then
          echo "  .$BASEDIR/scripts/run_split.sh -o $OUTNAME -p $THREADS $SPLIT_ARGS"
          bash $BASEDIR/scripts/run_split.sh -f $FQ0 -o $SRA_RUN -p $THREADS $SPLIT_ARGS & wait

        else
          echo "   ERROR: Neither paired or unpaired reads detected"
          # Update scheduler with fail message
          exit 1
        fi

        # Count output blocks
        N_paired=$( (ls *1.fq.* 2>/dev/null ) | wc -l)
        echo "    N paired-end fq-blocks: $N_paired"

        N_unpaired=$((ls *0.fq.* 2>/dev/null ) | wc -l)
        echo "    N unpaired   fq-blocks: $N_unpaired"
        echo ""

        # RUN UPLOAD ==============================================
        echo "  Running -- run_upload.sh --"
        echo "  ./scripts/run_upload.sh -k $S3_BUCKET -s $SRA_RUN $UL_ARGS"

        bash $BASEDIR/scripts/run_upload.sh \
             -k $S3_BUCKET \
             -s $SRA_RUN & wait

        echo "  Uploading complete."
        echo "  Status: DONE"

        # CLEAN-UP ================================================
        # rm may be unneccesary since container is shutdown
        cd $BASEDIR; rm -rf $WORKDIR/*

        # Update to scheduler -------------------------------------
        # ---------------------------------------------------------
        ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)
        echo "  $WORKERID - Job $ACC_ID is complete. Update scheduler."
        RESPONSE=$(curl -s "$SCHED/jobs/$ACC_ID?status=split_done&N_paired=$N_paired&N_unpaired=$N_unpaired")
    done
}


# SCRIPT CORE ===================================
echo "=========================================="
echo "                SERRATUS                  "
echo "=========================================="

# Default base directory is /home/serratus
cd $BASEDIR

boot_procedure
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
    INSTANCE_ACTION=$(curl -s $METADATA/spot/instance-action | jq -r .action)
    if [ "$INSTANCE_ACTION" == "terminate" ]; then
        echo "SPOT TERMINATION SIGNAL RECEIEVED."
        echo "Initiating shutdown procedures for all workers"

        for i in $(seq 1 "$WORKERS"); do
            kill -USR1 ${worker[i]} 2>/dev/null || true
        done
        break
    fi

    sleep 5
done

echo "============================"
echo "======= RUN COMPLETE ======="
echo "============================"
exit 0
