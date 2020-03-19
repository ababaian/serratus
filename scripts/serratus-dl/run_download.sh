#!/bin/bash
# run_download
#
set -e
# Script
SCRIPT_VERSION='0.1'
# Container: serratus-dl (0.1)
CONTAINER_VERSION='serratus-dl:0.1'

#Usage function
function usage {
  echo ""
  echo "Usage: run_download.sh -s <SRA Accession> [<D>] [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    SRA Accession for sratools"
  echo "    -s    SRA Accession for downloading SRA-archive file"
  echo ""
  echo "    Performance"
  echo "    -n CPU threads to use where applicable [1]"
  echo ""
  echo "    Arguments from serratus-dl"
  echo "    <D>   Arguments from running serratus-dl passed to this script (optional)"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_prefix> [SRA_ACCESSION]"
  echo ""
  echo "    Outputs: <output_prefix>.fq.0"
  echo "             <output_prefix>.fq.1"
  echo "             <output_prefix>.fq.2"
  echo ""
  echo "ex: ./run_download.sh -s SRR1337"
  exit 1
}

# PARSE INPUT =============================================
# SRA Accession -S
SRA=''
THREADS='1'

# Script Arguments -DPU
DL_ARGS=''

# Output options -do
WORKDIR="$PWD"
OUTNAME="$SRA"

while getopts s:n:D:d:o:h FLAG; do
  case $FLAG in
    # Input Options ---------
    s)
      SRA=$OPTARG
      ;;
    n)
      THREADS=$OPTARG
      ;;
    D)
      DL_ARGS=$OPTARG
      ;;
    # output options -------
    d)
      WORKDIR=$OPTARG
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

if [[ ( -z "$SRA" ) ]]
then
  echo "-s SRA_ACCESSION is required."
  usage
fi

# SCRIPT ==================================================

# Script --------------------
cd $WORKDIR

echo "      fastq-dump --split-e $SRA"
fastq-dump --split-e $SRA
if [ -f "$SRA.fastq" ]; then
    mv "$SRA.fastq" "$SRA"_0.fastq
fi

# TODO: Update to fasterq-dump
#echo "      fasterq-dump  $SRA -e $THREADS"
#fasterq-dump $SRA -e $THREADS

echo "      Download complete."
