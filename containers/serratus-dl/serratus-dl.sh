#!/bin/bash
set -eux
# ENTRYPOINT SCRIPT ===================
# serrtaus-dl
# =====================================
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
PIPE_VERSION="0.1.4"
AMI_VERSION='ami-0fdf24f2ce3c33243'
CONTAINER_VERSION='serratus-dl:0.1.4'

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
  echo "    Environment variables:"
  echo "      SCHEDULER: A url pointing to the scheduler.  This is overridden by -u"
  echo "      JOB_JSON: A JSON string containing the description of work to do"
  echo "ex: docker exec serratus-dl -u localhost:8000"
  exit 1
}
# Scheduler / Container Parameters

# SRA Accession -S
SRA=''
THREADS='1'

# AWS Options -ak
AWS_CONFIG='TRUE'
S3_BUCKET='serratus-public'

# Script Arguments -DPU
DL_ARGS=''
SPLIT_ARGS=''
UL_ARGS=''

# Output options -do
BASEDIR=${BASEDIR:-"/home/serratus"}
OUTNAME="$SRA"

while getopts u:w:n:s:ak:D:P:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    u)
      SCHEDULER=$OPTARG
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

if [ -z "$SCHEDULER" ]; then
    echo Please set SCHEDULER environment variable or use -k flag.
    exit 1
fi

# Parse SRA Accession ID
ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)

# Set up a error trap.  If something goes wrong unexpectedly, this will send
# a message to the scheduler before exiting.
function error {
    echo Error encountered.  Notifying the scheduler.
    curl -s -X POST "$SCHEDULER/jobs/split/$ACC_ID?state=split_err"
}
trap error ERR

SRA=$(echo $JOB_JSON | jq -r .sra_run_info.Run)
# Check inputs --------------
if [[ ( -z "$OUTNAME" ) ]];
    then OUTNAME="$SRA"
fi

# TODO: Allow the scheduler/main data-table to have arugments
# which will be passed on to the downloader scripts

# Run job --------------------------------------------------
# ----------------------------------------------------------
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )
WORKDIR=$BASEDIR/work/$RUNID
mkdir -p $WORKDIR; cd $WORKDIR

# Generate a unique downloader ID
DLID="$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )-\
$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 4 | head -n 1 )-\
$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 4 | head -n 1 )-\
$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 4 | head -n 1 )-\
$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 12 | head -n 1 )"
sed -i "s/52e8a8fe-0cac-4bf2-983a-3617cdba7df5/$DLID/g" /root/.ncbi/user-settings.mkfg

# Cleanup on exit, regardless of whether we encounter an error
function cleanup {
    rm -rf "$WORKDIR"
}
trap cleanup EXIT

echo "============================"
echo "    serratus-dl Pipeline    "
echo "============================"
echo " date:      $(date)"
echo " version:   $PIPE_VERSION"
echo " ami:       $AMI_VERSION"
echo " container: $CONTAINER_VERSION"
echo " worker-id: $WORKER_ID"
echo " run-id:    $RUNID"
echo " sra:       $SRA"
echo " S3 url:    $S3_BUCKET"
echo ""

# RUN DOWNLOAD ==============================================
echo "  Running -- run_dl-sra.sh --"
echo "  $BASEDIR/run_dl-sra.sh $SRA $DL_ARGS"

export BASEDIR
$BASEDIR/run_dl-sra.sh "$SRA" "$S3_BUCKET"

echo ''

# TODO Detecting paired / unpaired reads for the scheduler.
FQ1_PREFIX="fq-blocks/$SRA/$SRA.1.fq."
FQ2_PREFIX="fq-blocks/$SRA/$SRA.2.fq."
FQ3_PREFIX="fq-blocks/$SRA/$SRA.3.fq."

N_FQ1=$(aws s3api list-objects-v2 --bucket "$S3_BUCKET" --prefix "$FQ1_PREFIX" --query "length(Contents[])" || echo "0")
N_FQ2=$(aws s3api list-objects-v2 --bucket "$S3_BUCKET" --prefix "$FQ2_PREFIX" --query "length(Contents[])" || echo "0")
N_FQ3=$(aws s3api list-objects-v2 --bucket "$S3_BUCKET" --prefix "$FQ3_PREFIX" --query "length(Contents[])" || echo "0")

if [ "$N_FQ1" != "$N_FQ2" ]
then
  echo "Number of FQ1 files != FQ2 files.  Notifying scheduler"
  false # Handle error and quit
fi
  
if [[ ( "$N_FQ1" -gt 0 ) && ( "$N_FQ2" -gt 0 ) ]]
then
  paired_exists=true
  echo "  Paired-end reads detected"
else
  paired_exists=false
  echo "  Paired-end reads not-detected"
fi

if [[ ( "$N_FQ3" -gt 0 ) ]]
then
  unpaired_exists=true
  echo "  Unpaired reads detected"
else
  unpaired_exists=false
  echo "  unpaired reads not-detected"
fi

if [[ "$paired_exists" = false && "$unpaired_exists" = false ]]
then
  echo "Neither paired-end nor unpaired exist.  Notifying scheduler"
  false # Handle error and quit
fi

if [[ "$paired_exists" = true && "$unpaired_exists" = true ]]
then
  # Both paired-end and unpaired exist
  # use only paired-end data
  echo "   WARNING: Paired and Unpaired data detected"
  echo "            Using only Paired-End Reads"
  unpaired_exists=FALSE
  N_FQ3=0
fi

# Count output blocks
echo "    N paired-end fq-blocks: $N_FQ1"
echo "    N unpaired   fq-blocks: $N_FQ3"
echo ""

# Update to scheduler -------------------------------------
# ---------------------------------------------------------
ACC_ID=$(echo $JOB_JSON | jq -r .acc_id)
echo "  $WORKER_ID - Job $ACC_ID is complete. Update scheduler."
curl -s -X POST "$SCHEDULER/jobs/split/$ACC_ID?state=split_done&N_paired=$N_FQ1&N_unpaired=$N_FQ3"

echo "============================"
echo "======= RUN COMPLETE ======="
echo "============================"
exit 0
