#!/bin/bash
# run_resummarize
#
# Base: serratus-merge 

set -eu
PIPE_VERSION="0.3.0"

function usage {
  echo ""
  echo "Usage: run_resummarize.sh -s <SRR accession>' [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Parameters"
  echo "    -s    SRA Accession"
  echo "    -g    Genome identifier [cov3ma] (used for sumzer)"
  echo "    -L    S3 bucket path [s3://lovelywater2]"
  echo ""
  echo "    Merge Parameters"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo ""
  #echo "    Optional outputs"
  #echo "    -i    Flag. Generate bam.bai index file. Requires sort, otherwise false."
  #echo "    -f    Flag. Generate flagstat summary file"
  #echo "    -r    Flag. Sort final bam output (requires double disk usage)"
  #echo ""
  #echo ""
  echo "    Output options"
  echo "    -o    <output_filename_prefix> [Defaults to SRA_ACCESSION]"
  echo "    -O    Output S3 bucket [s3://serratus-bio]"
  echo ""
  echo "    Outputs a sorted Uploaded to s3: "
  echo "          <output_prefix>.bam, <output_prefix>.bam.bai, <output_prefix>.flagstat"
  echo ""
  echo "ex: bash run_resummarize.sh -s 'SRR123'"
  exit 1
}


# PARSE INPUT =============================================
# Generate random alpha-numeric for run-id
#RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Run Parameters
SRA=''
GENOME='cov3ma'
S3='s3://lovelywater2'

# Merge Options
THREADS='1'
#INDEX='negative'
#FLAGSTAT='negative'
#SORT='negative'

# Script Arguments -M
MERGE_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME=''
S3_OUT='s3://serratus-bio'


while getopts s:L:o:O:n:ifrh FLAG; do
  case $FLAG in
    s)
      SRA=$OPTARG
      ;;
    o)
      OUTNAME=$OPTARG
      ;;
    g)
      GENOME=$OPTARG
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
# Required parameters
if [ -z "$SRA" ]
then
  echo "-s <SRA / Outname> required"
  usage
fi

if [ -z "$OUTNAME" ]
then
  OUTNAME="$SRA"
fi

# Final Output Bam File name
OUTBAM="$OUTNAME.bam"


# SCRIPT ===================================================
# Command to run summarizer script
sumzer="$GENOME.sumzer.tsv"

if [ ! -e "$GENOME.sumzer.tsv" ]; then
        echo "  $GENOME.sumzer.tsv not found. Attempting download from"
        echo "  $S3/seq/$GENOME/$GENOME.sumzer.tsv"

        aws s3 cp $S3/seq/$GENOME/$GENOME.sumzer.tsv ./
fi


# Meta-data header for summary file
export SUMZER_COMMENT=$(echo sra="$SRA",genome="$GENOME",version=200818,date=$(date +%y%m%d-%R))

# Summary Comment / Meta-data
# usage: serratus_summarizer_flom.py InputSamFileName MetaTsvFilename SummaryFileName OutputSamFileName
#summarizer="python3 /home/serratus/serratus_summarizer.py /dev/stdin $sumzer $SRA.summary /dev/stdout"
summarizer="python3 /home/serratus/serratus_summarizer.py /dev/stdin $sumzer $SRA.summary /dev/null"

# Acquire + Run -----------------------
# Download bam file
aws s3 cp $S3/bam/$SRA.bam ./$SRA.unsorted.bam

# Summarize v2
samtools view $SRA.unsorted.bam | \
$summarizer 

# Sort
samtools sort -@ $THREADS $SRA.unsorted.bam >\
$OUTBAM

# index
samtools index $OUTBAM

# Upload ------------------------------
if [[ -s "$SRA.bam" ]]; then
  aws s3 cp --only-show-errors $SRA.bam $S3_OUT/bam/
fi

if [[ -s "$SRA.bam.bai" ]]; then
  aws s3 cp --only-show-errors $SRA.bam.bai $S3_OUT/bam/
fi

if [[ -s "$SRA.summary" ]]; then
  aws s3 cp --only-show-errors $SRA.summary $S3_OUT/summary/
fi

# Clean-up
rm $SRA.unsorted.bam $OUTBAM $OUTBAM.bai $SRA.summary

# end of script