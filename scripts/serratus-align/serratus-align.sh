#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-align
# =====================================
set -e
#
# Base: serratus-aligner (0.1)
# Amazon Linux 2 with Docker
# AMI : aami-0fdf24f2ce3c33243
# login: ec2-user@<ipv4>
# base: 9 Gb
# Container: serratus-align:0.1
#

# Test cmd: ./serratus-align.sh 
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
  echo "Usage: docker exec serratus-align -u [localhost:8000] -k <s3://bucket/bam-blocks> [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Fields"
  echo "    -u    <IP>:<PORT> to query for scheduler webserver [localhost:8000]"
  echo "    -k    S3 bucket path for bam-blocks [s3://serratus-public/bam-blocks]"
  echo ""
  echo "    Alignment Parameters"
  echo "    -a    Alignment software [bowtie2] (no alt. currently)"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo ""
  echo "    Manual overwrites"
  echo "    -s    (N/A) SRA Accession for direct-to-pipeline access (from scheduler)"
  echo "    -q    (N/A) fastq-block(s) (s3 URL) to run aligner on"
  echo "    -g    (N/A) Genome Reference name (from scheduler)"
  echo "           genome index at s3://serratus-public/resources/<GENOME>/*"
  echo ""
  echo "    AWS / S3 Bucket parameters"
  echo "    -w    Flag. Imports AWS IAM from this ec2-instance for container"
  echo "          EC2 instance must have been launched with correct IAM"
  echo "          (No alternative yet, hard set to TRUE)"
  echo ""
  echo "    Arguments as string to pass to downstream scripts"
  echo "    -A    String of arguments to pass to 'run_bowtie2.sh <ARG>'"
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
S3_BUCKET='s3://serratus-public/bam-blocks'

# Aligner
ALIGNER='bowtie2'
THREADS='1'

# Inputs (manual overwrite)
SRA=''
GENOME=''
S3_FQ0=''
S3_FQ1=''
S3_FQ2=''

# AWS Options -ak
AWS_CONFIG='TRUE'

# Script Arguments -AU
ALIGN_ARGS=''
UL_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts u:k:a:n:s:g:0:1:2:A:U:d:owh FLAG; do
  case $FLAG in
    # Scheduler Options ---------
    u)
      SCHED=$OPTARG
      ;;
    k)
      S3_BUCKET=$OPTARG
      ;;
    # Aligner Options ---------
    a)
      #ALIGNER=$OPTARG
      ALIGNER="bowtie2"
      ;;
    n)
      THREADS=$OPTARG
      ;;
    s)
      SRA=$OPTARG
      ;;
    g)
      GENOME=$OPTARG
      ;;
    0)
      S3_FQ0=$OPTARG
      ;;
    1)
      S3_FQ1=$OPTARG
      ;;
    2)
      S3_FQ2=$OPTARG
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
#if [[ ( -z "$OUTNAME" ) ]]; then OUTNAME="$SRA"; fi

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
# SRA: SRA accession / run-name
# PAIRED: [true / false] is there paired-end fq data to process, else un-paired
# FQ0: S3 path to unpaired fq-block
# FQ1: S3 path to paired fq-block 1/2
# FQ2: S3 path to paired fq-block 2/2
# BL_N: Block number
# GENOME: genome reference name
# ALIGN_ARG: Arguments to pass to aligner
# 
## Read Group meta-data (required for GATK) -LISPF
# RGLB: Library Name / SRA Accession
# RGID: Run ID / SRA Accession (?)
# RGSM: Sample ID / BioSample
# RGPO: Population / Experiment
# RGPL: Platform [ILLUMINA]

  while true; do
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
      *)        echo "  $WORKERID - ERROR: Unknown State received."
        exit 1
        echo "  $WORKERID - ERROR: Unknown State received."
        exit 1
    esac

    # Parse Job JSON for run parameters
    SRA=$(jq    -r .acc_id     <(echo $JOB_JSON))

    PAIRED=$(jq -r .is_paired  <(echo $JOB_JSON))
    BL_N=$(jq   -r .block_n    <(echo $JOB_JSON))
    S3_FQ0=$(jq -r .s3_fq0     <(echo $JOB_JSON))
    S3_FQ1=$(jq -r .s3_fq1     <(echo $JOB_JSON))
    S3_FQ2=$(jq -r .s3_fq2     <(echo $JOB_JSON))

    ALIGN_ARGS=$(jq -r .align_args <(echo $JOB_JSON))
    GENOME=$(jq     -r .genome     <(echo $JOB_JSON))

    RGLB=$(jq -r .acc_id        <(echo $JOB_JSON))
    RGID=$(jq -r .acc_id        <(echo $JOB_JSON))
    RGSM=$(jq -r .biosample_id  <(echo $JOB_JSON))
    RGPO=$(jq -r .experiment_id <(echo $JOB_JSON))
    RGPL=$(jq -r .platform      <(echo $JOB_JSON))

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
    GENDIR=$BASEDIR/$GENOME

    echo "============================"
    echo "  serratus-align Pipeline   "
    echo "============================"
    echo " date:      $(date)"
    echo " version:   $PIPE_VERSION"
    echo " ami:       $AMI_VERSION"
    echo " container: $CONTAINER_VERSION"
    echo " run-id:    $RUNID"
    echo " sra:       $SRA"
    echo " block:     $BL_N"
    echo " genome:    $GENOME"
    echo " align arg: $ALIGN_ARGS"
    echo " rg:        --rglb $RGLB --rgid $RGID --rgsm $RGSM --rgpo $RGPO --rgpl $RGPL"


# DOWNlOAD GENOME =========================================
    if [ -d $GENDIR ]
    then
      echo "  $GENDIR found."
    else
      echo "  $GENDIR not found. Attempting download from"
      echo "  s3://serratus-public/resources/$GENOME"
      mkdir -p $GENEDIR; cd $GENEDIR

      aws s3 cp --recursive s3://serratus-public/$GENOME/ $GENEDIR/

      if [[ -e "$GENOME.fa" && -e "$GENOME.1.bt2" ]]
      then
        echo " ERROR: $GENOME.fa or $GENOME.1.bt2 index not found"
        exit 1
      else
        echo "  genome download complete"
      fi
    fi

    # Link genome files to workdir
    cd $WORKDIR
    ln -s $GENEDIR/* ./

# DOWNlOAD FQ Files =======================================

    if [[ "$PAIRED" = true ]]
    then
      echo "  Downloading Paired-end fq-block data..."
      aws s3 cp $S3_FQ1 ./
      aws s3 cp $S3_FQ2 ./

      FQ1=$(basename $S3_FQ1)
      FQ2=$(basename $S3_FQ2)

    else
      echo "  Downloading unpaired fq-block data..."
      aws s3 cp $S3_FQ0 ./

      FQ0=$(basename $S3_FQ0)
    fi

    ## Test data
    # S3_FQ1='s3://serratus-public/fq-blocks/SRR11166696/SRR11166696.1.fq.0000000000.gz'
    # S3_FQ2='s3://serratus-public/fq-blocks/SRR11166696/SRR11166696.2.fq.0000000000.gz'
    # RGLB='tmp'; RGID='tmp2'; RGSM='tmp3'; RGPO='tmp4'

# RUN ALIGN ===============================================
    echo "  Running -- run_bowtie2.sh --"

    if [[ "$PAIRED" = true ]]
    then
      echo "  bash $BASEDIR/scripts/run_bowtie2.sh " &&\
      echo "    -1 $FQ1 -2 $FQ2 -x $GENOME" &&\
      echo "    -o $SRA.$BL_N -p $THREADS $ALIGN_ARGS" &&\
      echo "    -L $RGLB -I $RGID -S $RGSM -P $RGPO"

      bash $BASEDIR/scripts/run_bowtie2.sh \
        -1 $FQ1 -2 $FQ2 -x $GENOME \
        -o $SRA.$BL_N -p $THREADS $ALIGN_ARGS \
        -L $RGLB -I $RGID -S $RGSM -P $RGPO & wait
    else
      echo "  bash $BASEDIR/scripts/run_bowtie2.sh " &&\
      echo "    -0 $FQ0 -x $GENOME" &&\
      echo "    -o $SRA.$BL_N -p $THREADS $ALIGN_ARGS" &&\
      echo "    -L $RGLB -I $RGID -S $RGSM -P $RGPO"

      bash $BASEDIR/scripts/run_bowtie2.sh \
        -0 $FQ0 -x $GENOME \
        -o $SRA.$BL_N -p $THREADS $ALIGN_ARGS \
        -L $RGLB -I $RGID -S $RGSM -P $RGPO & wait
    fi

# RUN UPLOAD ==============================================
    echo "  Uploading bam-block data..."
    echo "  $SRA.$BL_N.bam"

    aws s3 cp $SRA.$BL_N.bam $S3_BUCKET/$SRA/

    echo "  Status: DONE"

# CLEAN-UP ================================================
  # Update to scheduler -------------------------------------
  # ---------------------------------------------------------
    # Tell the scheduler we're done
    echo "  $WORKERID - Job $SRA is complete. Update scheduler."
    RESPONSE=$(curl -s "$SCHED/finish_split_job?job_id=$SRA&status=align_done&N=$N")

    cd $BASEDIR; rm -rf $WORKDIR/*

  # Free up fq-blocks from s3
  if [[ "$PAIRED" = true ]]
    then
      aws s3 rm $S3_FQ1
      aws s3 rm $S3_FQ2
    else
      aws s3 rm $S3_FQ0
    fi
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
