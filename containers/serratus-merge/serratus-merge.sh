#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-merge
# =====================================
set -eu
#
# Amazon Linux 2 with Docker
# AMI : aami-0fdf24f2ce3c33243
# login: ec2-user@<ipv4>
# base: 50 Gb
PIPE_VERSION="0.3.0"
CONTAINER_VERSION='serratus-merge:0.3.0'

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
  echo "    -i    Flag. Generate bam.bai index file"
  echo "    -f    Flag. Generate flagstat summary file"
  echo "    -r    Flag. Sort final bam output (requires double disk usage)"
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
  exit 0
}

# PARSE INPUT =============================================
# Scheduler / Container Parametersterraform apply -auto-approve
S3_BUCKET=${S3_BUCKET:-'serratus-public'}

# Merge Options
THREADS='1'
INDEX='false'
FLAGSTAT='false'
SORT='false'

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

while getopts u:k:n:s:b:g:M:d:o:whifr FLAG; do
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
      INDEX="true"
      ;;
    f)
      FLAGSTAT="true"
      ;;
    r)
      SORT="true"
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
    false
    exit 0
fi

# SCRIPT CORE ===================================

# Default base directory is /home/serratus
cd $BASEDIR

# Query for job --------------------------------------------
# ----------------------------------------------------------
# DATA NEEDED FROM SCHEDULER
# SRA   :   SRA accession / run-name / output name
# EXP   :   Experiment ID
# GENOME:   genome reference name
# S3_BAM:   S3 download directory to bam-blocks
# BLOCKS:   Number of bam-blocks expected

ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)

# Set up a error trap.  If something goes wrong unexpectedly, this will send
# a message to the scheduler before exiting.
function error {
    echo Error encountered.  Notifying the scheduler.
    curl -s -X POST "$SCHEDULER/jobs/merge/$ACC_ID?state=merge_err" > /dev/null
    exit 0 # Already told the calling script.
}
trap error ERR

# Parse Job JSON for run parameters
SRA=$(echo $JOB_JSON    | jq -r .sra_run_info.Run)
PAIRED=$(echo $JOB_JSON | jq -r .contains_paired)
BLOCKS=$(echo $JOB_JSON | jq -r .blocks)
GENOME=$(echo $JOB_JSON | jq -r .genome)

#TODO IMPLEMENT EXP NAME FROM SCHEDULER
#EXP=$(echo $JOB_JSON     | jq -r .expid)

MERGE_ARGS=$(echo $JOB_JSON | jq -r .merge_args)

# Run job --------------------------------------------------
# ----------------------------------------------------------
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )
WORKDIR=$BASEDIR/work/$RUNID
mkdir -p $WORKDIR; cd $WORKDIR

# Download pan-genome summary data
if [ ! -f "$BASEDIR/$GENOME.sumzer.tsv" ]; then
  aws s3 cp --only-show-errors "s3://serratus-public/seq/$GENOME/$GENOME.sumzer.tsv" $BASEDIR/
fi

cp $BASEDIR/"$GENOME".sumzer.tsv $WORKDIR/ # copy summarizer tables into runid

S3_BAM=s3://$S3_BUCKET/bam-blocks/$SRA
S3_FQ=s3://$S3_BUCKET/fq-blocks/$SRA

# COUNT + DL BLOCKS =======================================
# Query S3_OUT and count number of matching files
# which should match BLOCKS parameter from scheduler

BL_COUNT=$(aws s3 ls $S3_BAM/ | wc -l | cut -f1 -d' ' -)

if [[ "$BLOCKS" != "$BL_COUNT" ]]
 then
   echo "  ERROR: Number of bam-blocks in $S3_OUT: $BL_COUNT"
   echo "         is not equal to the expected $BLOCKS"
   echo ""
   echo "  aws s3 ls $S3_BAM/"
   echo ""
   aws s3 ls $S3_BAM/
   false
   exit 1
 fi

# Download bam-blocks
aws s3 cp --only-show-errors --recursive $S3_BAM ./

# RUN MERGE ===============================================
bash $BASEDIR/run_merge.sh -s $SRA -b "*bam" -x $GENOME


# RUN UPLOAD ==============================================
if [[ -s "$SRA.bam" ]]; then
  aws s3 cp --only-show-errors $SRA.bam $S3_OUT/bam/
fi

# Not generated by default
if [[ -s "$SRA.bam.bai" ]]; then
  aws s3 cp --only-show-errors $SRA.bam.bai $S3_OUT/bai/
fi

if [[ -s "$SRA.flagstat" ]]; then
  aws s3 cp --only-show-errors $SRA.flagstat $S3_OUT/flagstat/
fi    

if [[ -s "$SRA.summary" ]]; then
  aws s3 cp --only-show-errors $SRA.summary $S3_OUT/summary/
fi

if [[ -s "$SRA.th" ]]; then
  aws s3 cp --only-show-errors $SRA.th $S3_OUT/th/
fi

# CLEAN-UP ================================================
# Tell the scheduler we're done
curl -X POST -s "$SCHEDULER/jobs/merge/$ACC_ID?state=merge_done" > /dev/null

cd $BASEDIR; rm -rf $WORKDIR/*

# Free up fq-blocks and bam-blocks from s3
aws s3 rm --only-show-errors --recursive $S3_BAM
aws s3 rm --only-show-errors --recursive $S3_FQ

exit 0

