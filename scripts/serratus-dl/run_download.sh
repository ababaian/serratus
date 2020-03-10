#!/bin/bash
# run_download
#
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
  echo "    Arguments from serratus-dl"
  echo "    <D>   Arguments from running serratus-dl passed to this script (optional)"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [/home/serratus]"
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

# Script Arguments -DPU
DL_ARGS=''

# Output options -do
WORKDIR="/home/serratus"
OUTNAME="$SRA"

while getopts s:ak:D:P:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    s)
      SRA=$(readlink -f $OPTARG)
      ;;
    D)
      DL_ARGS=$(readlink -f $OPTARG)
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

# Script --------------------
cd $WORKDIR

echo "      fastq-dump --split-e $SRA"
fastq-dump --split-e $SRA