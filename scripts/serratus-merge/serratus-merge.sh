#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-merge
# =====================================
set -e
#
# Base: serratus-merger (0.1)
# Amazon Linux 2 with Docker
# AMI : aami-0fdf24f2ce3c33243
# login: ec2-user@<ipv4>
# base: 9 Gb
# Container: serratus-merge:0.1
#m

# Test cmd: ./serratus-merge.sh 
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
CONTAINER_VERSION='serratus-merge:0.1'

# Usage
function usage {
  echo ""
  echo "Usage: docker exec serratus-merge -u [localhost:8000] -k <s3://bucket/out/path> [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Fields"
  echo "    -u    <IP>:<PORT> to query for scheduler webserver [localhost:8000]"
  echo "    -k    S3 bucket path for /bam, /bai, /flagstat [s3://serratus-public/out]"
  echo ""
  echo "    Merge Parameters"
  echo "    -i    Flag. Do not generate bam.bai index file"
  echo "    -f    Flag. Do not generate flagstat summary file"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo ""
  echo "    Manual overwrites"
  echo "    -s    (N/A) SRA Accession for direct-to-pipeline access (from scheduler)"
  echo "    -b    (N/A) bam-block(s) (s3 URL) to run merger on"
  echo "    -g    (N/A) Genome Reference name (from scheduler)"
  echo "           genome index at s3://serratus-public/resources/<GENOME>/*"
  echo ""
  echo "    AWS / S3 Bucket parameters"
  echo "    -w    Flag. Imports AWS IAM from this ec2-instance for container"
  echo "          EC2 instance must have been launched with correct IAM"
  echo "          (No alternative yet, hard set to TRUE)"
  echo ""
  echo "    Arguments as string to pass to downstream scripts"
  echo "    -M    String of arguments to pass to 'run_merge.sh <ARG>'"
  echo ""
  echo "    Output options"
  echo "    -d    Base directory in container [/home/serratus/]"
  echo "    -o    <output_filename_prefix> [Defaults to SRA_ACCESSION]"
  echo ""
  echo "    Outputs Uploaded to s3: "
  echo "          <output_prefix>.bam ... <output_prefix>.fq.yyyy"
  echo ""
  echo "ex: docker exec serratus-dl -u localhost:8000"
  exit 1
}

# PARSE INPUT =============================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Scheduler / Container Parameters
SCHED='localhost:8000'
S3_BUCKET='s3://serratus-public/out'

# Merge Options
THREADS='1'
INDEX='true'
FLAGSTAT='true'

# Inputs (manual overwrite)
SRA=''
S3_OUT=''
GENOME=''

# AWS Options -ak
AWS_CONFIG='TRUE'

# Script Arguments -AU
MERGE_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts u:k:n:s:b:g:M:d:o:whif FLAG; do
  case $FLAG in
    # Scheduler Options -----
    u)
      SCHED=$OPTARG
      ;;
    k)
      S3_BUCKET=$OPTARG
      ;;
    # Merge Options ---------
    n)
      THREADS=$OPTARG
      ;;
    i)
      INDEX="false"
      ;;
    f)
      FLAGSTAT="false"
      ;;
    # Manual Overwrite ------
    s)
      SRA=$OPTARG
      ;;
    g)
      GENOME=$OPTARG
      ;;
    b)
      S3_OUT=$OPTARG
      ;;
    w)
      AWS_CONFIG='TRUE'
      ;;
    # SCRIPT ARGUMENTS -------
    A)
      ALIGN_ARGS=$OPTARG
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
if [ -z $S3_BAM ]
then
  echo "(-o) output prefix is required."
  usage
fi

# FUNCTIONS ===============================================
# Worker process.  Run in a loop, grabbing data from the scheduler and
# processing it.

function boot_procedure {
  # AUTHENTICATE AWS ========================================
  echo "  Authenticating AWS credentials"

  # Create aws key (read from ec2 IAM) (run at start-up of container)
  export AWS_ACCESSKEYID=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.AccessKeyId' )
  export AWS_SECRETKEY=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.SecretAccessKey' )
  
  # Create a key file
  export AWS_KEY="$WORKDIR/key.csv"
  echo User name,Password,Access key ID,Secret access key,Console login link > $AWS_KEY
  echo default,,$AWS_ACCESSKEYID,$AWS_SECRETKEY, >> $AWS_KEY
  chmod 400 $AWS_KEY

  # Pass credentials to sra-toolkit via vdb-config
  # Current version requires a manual
  # vdb-config -i initialization
  # ugly hack is to copy a blank config file and bypass this
  # TODO: Move vdb-config -i hack to Dockerfile/installation
  mkdir -p /root/.ncbi
  aws s3 cp s3://serratus-public/VDB_user-settings.mkfg /root/.ncbi/user-settings.mkfg
  chmod 500 /root/.ncbi/user-settings.mkfg
  vdb-config --accept-aws-charges yes \
    --report-cloud-identity yes \
    --set-aws-credentials $AWS_KEY

  # Download AWS S3 test token
  aws s3 cp s3://serratus-public/aws-test-token.jpg $WORKDIR/aws-test-token.jpg

  if [ ! -s "$AWS_KEY" ]
  then
    echo "    ERROR: AWS Test did not download"
    echo "    Ensure the EC2 instance has correct IAM permissions"
    echo "      - requires S3 Read/Write"
    usage
  fi

  echo '  ...authentication token was download successfully'
  echo ""
}

function wait_for_scheduler {
    while true; do
        curl -s "$SCHED/status" && break
        sleep 5
    done
}

function terminate_handler {
    echo "    $SRA was terminated without completing. Reset status."
    echo "    In trap $(date -In)"
    # Tell server to reset this job to a "new" state, since we couldn't
    # finish processing it.
    curl -s "$SCHED/finish_split_job?job_id=$SRA&status=new"
}

function main_loop {
  trap terminate_handler SIGUSR1
  # Note: The "& wait"s are important.  Without them, bash will wait for
  # the command to finish before executing its traps.  When we use "& wait",
  # the command will recieve the same trap (killing it), and then run our
  # trap handler, which tells the server our job failed.
  # what if we & wait on the whole main_loop function instead? Test this.
  WORKER_ID=$1

# Query for job --------------------------------------------
# ----------------------------------------------------------
# DATA NEEDED FROM SCHEDULER
# SRA:    SRA accession / run-name / output name
# S3_BAM: S3 download directory to bam-blocks
# BL_N:   Number of bam-blocks expected

  while true; do
    echo "$WORKERID - Requesting merge job from Scheduler..."
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
      *)        echo "  $WORKERID - ERROR: Unknown State received."
        exit 1
        echo "  $WORKERID - ERROR: Unknown State received."
        exit 1
    esac

    # Parse Job JSON for run parameters
    SRA=$(jq    -r .acc_id     <(echo $JOB_JSON))
    S3_BAM=$(jq -r .s3_bam     <(echo $JOB_JSON))
    BL_N=$(jq   -r .block_n    <(echo $JOB_JSON))
    
    MERGE_ARGS=$(jq -r .merge_args <(echo $JOB_JSON))
    
    # TODO: Allow the scheduler/main data-table to have arugments
    # which will be passed on to the downloader scripts
    
    # TODO: Sort out how multiple 'workers' per instance
    # can be used for optimal merging.

# Run job --------------------------------------------------
# ----------------------------------------------------------
    # Generate random alpha-numeric for run-id
    RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )
    WORKDIR=$BASEDIR/$RUNID
    mkdir -p $WORKDIR; cd $WORKDIR
    GENDIR=$BASEDIR/$GENOME

    echo "============================"
    echo "  serratus-merge Pipeline   "
    echo "============================"
    echo " date:      $(date)"
    echo " version:   $PIPE_VERSION"
    echo " ami:       $AMI_VERSION"
    echo " container: $CONTAINER_VERSION"
    echo " run-id:    $RUNID"
    echo " sra:       $SRA"
    echo " bam path:  $S3_BAM"
    echo " blocks:    $BL_N"
    echo " make bai:  $INDEX"
    echo " make flagstat: $FLAGSTAT"
    echo " merge arg: $MERGE_ARGS"

# COUNT + DL BLOCKS =======================================
    # Query S3_OUT and count number of matching files
    # which should match BL_N parameter from scheduler

    # If there is a terminal slash on S3_BAM remove it
    S3_BAM=$(echo $S3_BAM | sed -i 's/\/$//g' -)
    BL_COUNT=$(aws s3 ls $S3_BAM/ | wc -l | cut -f1 -d' ' -)

    if [[ "$BL_N" != "BL_COUNT" ]]
     then
       echo "  ERROR: Number of bam-blocks in $S3_OUT"
       echo "         is not equal to the expected $BL_N"
       echo ""
       echo "  aws s3 ls $S3_BAM/"
       echo ""
       aws s3 ls $S3_BAM/
       exit 1
     fi

     # Download bam-blocks
     aws s3 cp $S3_BAM ./

# RUN MERGE ===============================================
    echo "  Running -- run_merge.sh --"
    echo ""
    bash $BASEDIR/scripts/run_merge.sh -s $SRA -b "*bam" 


# RUN UPLOAD ==============================================
    echo "  Uploading final bam data..."

    if [[ -s "$SRA.bam" ]]
    then
      echo "  ...$SRA.bam detected. Uploading."
      aws s3 cp $SRA.bam $S3_OUT/bam/
    fi
    
    if [[ -s "$SRA.bam.bai" ]]
    then
      echo "  ...$SRA.bam.bai detected. Uploading."
      aws s3 cp $SRA.bam.bai $S3_OUT/bam/
    fi

    if [[ -s "$SRA.flagstat" ]]
    then
      echo "  ...$SRA.flagstat detected. Uploading."
      aws s3 cp $SRA.flagstat $S3_OUT/flagstat/
    fi    
    echo "  Status: DONE"

# CLEAN-UP ================================================
  # Update to scheduler -------------------------------------
  # ---------------------------------------------------------
    # Tell the scheduler we're done
    echo "  $WORKERID - Job $SRA is complete. Update scheduler."
    RESPONSE=$(curl -s "$SCHED/finish_split_job?job_id=$SRA&status=merge_done&N=$N")

    cd $BASEDIR; rm -rf $WORKDIR/*

    # Free up bam-blocks from s3
    aws s3 rm --recursive $S3_BAM
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

echo "============================"
echo "======= RUN COMPLETE ======="
echo "============================"
exit 0


# Spot Operations ===============================
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
