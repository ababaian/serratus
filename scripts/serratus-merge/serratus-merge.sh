#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-merge
# =====================================
set -eu
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
S3_BUCKET=${S3_BUCKET:-'serratus-public'}

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
      SCHEDULER=$OPTARG
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

if [ -z "$SCHEDULER" ]; then
    echo Please set SCHEDULER environment variable or use -k flag.
    exit 1
fi

# SCRIPT CORE ===================================
echo "=========================================="
echo "                SERRATUS                  "
echo "=========================================="

# Default base directory is /home/serratus
cd $BASEDIR

# Query for job --------------------------------------------
# ----------------------------------------------------------
# DATA NEEDED FROM SCHEDULER
# SRA:    SRA accession / run-name / output name
# S3_BAM: S3 download directory to bam-blocks
# BLOCKS:   Number of bam-blocks expected

ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)

# Set up a error trap.  If something goes wrong unexpectedly, this will send
# a message to the scheduler before exiting.
function error {
    echo Error encountered.  Notifying the scheduler.
    curl -s -X POST "$SCHEDULER/jobs/merge/$ACC_ID?state=merge_err"
    exit 0 # Already told the calling script.
}
trap error ERR

# Parse Job JSON for run parameters
SRA=$(echo $JOB_JSON    | jq -r .sra_run_info.Run)
PAIRED=$(echo $JOB_JSON | jq -r .contains_paired)
BLOCKS=$(echo $JOB_JSON | jq -r .blocks)

MERGE_ARGS=$(echo $JOB_JSON | jq -r .merge_args)

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

S3_BAM=s3://$S3_BUCKET/bam-blocks/$SRA

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
echo " blocks:    $BLOCKS"
echo " make bai:  $INDEX"
echo " make flagstat: $FLAGSTAT"
echo " merge arg: $MERGE_ARGS"

# COUNT + DL BLOCKS =======================================
# Query S3_OUT and count number of matching files
# which should match BLOCKS parameter from scheduler

BL_COUNT=$(aws s3 ls $S3_BAM/ | wc -l | cut -f1 -d' ' -)

if [[ "$BLOCKS" != "BL_COUNT" ]]
 then
   echo "  ERROR: Number of bam-blocks in $S3_OUT"
   echo "         is not equal to the expected $BLOCKS"
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
echo "  $WORKER_ID - Job $SRA is complete. Update scheduler."
curl -X POST -s "$SCHEDULER/jobs/merge/$ACC_ID?state=merge_done"

cd $BASEDIR; rm -rf $WORKDIR/*

# Free up bam-blocks from s3
aws s3 rm --recursive $S3_BAM

echo "============================"
echo "======= RUN COMPLETE ======="
echo "============================"
exit 0

