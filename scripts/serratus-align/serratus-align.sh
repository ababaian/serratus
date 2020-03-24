#!/bin/bash
# ENTRYPOINT SCRIPT ===================
# serrtaus-align
# =====================================
set -eux
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
CONTAINER_VERSION='serratus-align:0.1'

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
S3_BUCKET=${S3_BUCKET:-'serratus-public'}

# Aligner
ALIGNER='bowtie2'
THREADS='1'

# Inputs (manual overwrite)
SRA=''
GENOME=''
S3_FQ1=''
S3_FQ2=''
S3_FQ3=''

# AWS Options -ak
AWS_CONFIG='TRUE'

# Script Arguments -AU
ALIGN_ARGS=''
UL_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts u:k:a:n:s:g:0:1:2:A:U:d:o:wh FLAG; do
  case $FLAG in
    # Scheduler Options ---------
    u)
      SCHEDULER=$OPTARG
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
    1)
      S3_FQ1=$OPTARG
      ;;
    2)
      S3_FQ2=$OPTARG
      ;;
    3)
      S3_FQ3=$OPTARG
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

echo "=========================================="
echo "                SERRATUS                  "
echo "=========================================="
cd $BASEDIR

# Query for job --------------------------------------------
# ----------------------------------------------------------
# DATA NEEDED FROM SCHEDULER
# SRA: SRA accession / run-name
# PAIRED: [true / false] is there paired-end fq data to process, else un-paired
# FQ1: S3 path to paired fq-block 1/2
# FQ2: S3 path to paired fq-block 2/2
# FQ3: S3 path to unpaired fq-block
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

echo "$WORKER_ID - Requesting align job from Scheduler..."

# Parse Job JSON for run parameters
BLOCK_ID=$(echo $JOB_JSON | jq -r .block_id)

# Set up a error trap.  If something goes wrong unexpectedly, this will send
# a message to the scheduler before exiting.
function error {
    curl -s -X POST "$SCHEDULER/jobs/align/$BLOCK_ID?state=fail"
    exit 0 # Already told the calling script.
}
trap error ERR

SRA=$(echo $JOB_JSON      | jq -r .sra_run_info.Run)
PAIRED=$(echo $JOB_JSON   | jq -r .contains_paired)
BL_N=$(echo $JOB_JSON     | jq -r .n)

ALIGN_ARGS=$(echo $JOB_JSON | jq -r .align_args)
GENOME=$(echo $JOB_JSON     | jq -r .genome)

RGLB="$SRA"
RGID="$SRA"
RGSM=$(echo $JOB_JSON | jq -r .sra_run_info.BioSample)
RGPO=$(echo $JOB_JSON | jq -r .sra_run_info.Experiment)
RGPL=$(echo $JOB_JSON | jq -r .sra_run_info.Platform)

S3_FQ1=$(printf 's3://%s/fq-blocks/%s/%s.1.fq.%010d.gz' "$S3_BUCKET" "$SRA" "$SRA" "$BL_N")
S3_FQ2=$(printf 's3://%s/fq-blocks/%s/%s.2.fq.%010d.gz' "$S3_BUCKET" "$SRA" "$SRA" "$BL_N")
S3_FQ3=$(printf 's3://%s/fq-blocks/%s/%s.3.fq.%010d.gz' "$S3_BUCKET" "$SRA" "$SRA" "$BL_N")

# TODO: Allow the scheduler/main data-table to have arugments
# which will be passed on to the aligner scripts

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
# This is not threadsafe!  For now let's just put it in a critical section.
# Using a recipe from "man flock" which appears to work.
(
    flock 200
    if [ -d $GENDIR ]; then
        echo "  $GENDIR found."
    else
        echo "  $GENDIR not found. Attempting download from"
        echo "  s3://serratus-public/resources/$GENOME"
        mkdir -p $GENDIR; cd $GENDIR

        aws s3 cp --recursive s3://serratus-public/$GENOME/ $GENDIR/

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
    ln -s $GENDIR/* ./
) 200> .genome-lock

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
  aws s3 cp $S3_FQ3 ./

  FQ3=$(basename $S3_FQ3)
fi

## Test data
# S3_FQ1='s3://serratus-public/fq-blocks/SRR11166696/SRR11166696.1.fq.0000000000.gz'
# S3_FQ2='s3://serratus-public/fq-blocks/SRR11166696/SRR11166696.2.fq.0000000000.gz'
# RGLB='tmp'; RGID='tmp2'; RGSM='tmp3'; RGPO='tmp4'

# RUN ALIGN ===============================================
echo "  Running -- run_bowtie2.sh --"

if [[ "$PAIRED" = true ]]
then
  echo "  bash $BASEDIR/run_bowtie2.sh " &&\
  echo "    -1 $FQ1 -2 $FQ2 -x $GENOME" &&\
  echo "    -o $SRA.$BL_N -p $THREADS -a $ALIGN_ARGS" &&\
  echo "    -L $RGLB -I $RGID -S $RGSM -P $RGPO"

  bash $BASEDIR/run_bowtie2.sh \
    -1 $FQ1 -2 $FQ2 -x $GENOME \
    -o $SRA.$BL_N -p $THREADS -a $ALIGN_ARGS \
    -L $RGLB -I $RGID -S $RGSM -P $RGPO
else
  echo "  bash $BASEDIR/run_bowtie2.sh " &&\
  echo "    -0 $FQ3 -x $GENOME" &&\
  echo "    -o $SRA.$BL_N -p $THREADS -a $ALIGN_ARGS" &&\
  echo "    -L $RGLB -I $RGID -S $RGSM -P $RGPO"

  bash $BASEDIR/run_bowtie2.sh \
    -0 $FQ3 -x $GENOME \
    -o $SRA.$BL_N -p $THREADS -a $ALIGN_ARGS \
    -L $RGLB -I $RGID -S $RGSM -P $RGPO
fi

# RUN UPLOAD ==============================================
echo "  Uploading bam-block data..."
echo "  $SRA.$BL_N.bam"

aws s3 cp $SRA.$BL_N.bam s3://$S3_BUCKET/bam-blocks/$SRA/

echo "  Status: DONE"

# CLEAN-UP ================================================
# Update to scheduler -------------------------------------
# ---------------------------------------------------------
# Tell the scheduler we're done
echo "  $WORKER_ID - Job $SRA is complete. Update scheduler."
curl -X POST -s "$SCHEDULER/jobs/align/$BLOCK_ID&state=done&N=$BL_N"

cd $BASEDIR; rm -rf $WORKDIR/*

# Free up fq-blocks from s3
if [[ "$PAIRED" = true ]]; then
    aws s3 rm $S3_FQ1
    aws s3 rm $S3_FQ2
else
    aws s3 rm $S3_FQ3
fi

echo "============================"
echo "======= RUN COMPLETE ======="
echo "============================"
exit 0
