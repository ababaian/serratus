#!/bin/bash
# run_upload
#
# Script
SCRIPT_VERSION='0.1'
# Container: serratus-dl (0.1)
CONTAINER_VERSION='serratus-dl:0.1'

#Usage function
function usage {
  echo ""
  echo "Usage: run_upload.sh -k <S3_BUCKET> -s <SRA/Name> [<-U>] [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    S3 Bucket parameters"
  echo "    -k    S3 URL to upload data to [s3://serratus-public/fq-blocks]"
  echo "    -s    SRA/Accession name"
  echo "    -r    REGEX string to match files for uploading [\"*.[123].fq*\"]"
  echo ""
  echo "    Arguments from serratus-dl"
  echo "    <-U>   Arguments as string passed from serratus-dl to"
  echo "           this script (optional)"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [PWD]"
  echo "    -o    <output_prefix> [-s option]"
  echo ""
  echo "ex: ./run_upload.sh -k s3://serratus-public/fq-blocks -s SRA1337"
  echo ""
  echo "    Will upload all files matching *.[123].fq* in <working dir> to"
  echo "    s3://serratus-public/fq-blocks/SRA1337/"
  echo ""
  exit 1
}

# PARSE INPUT =============================================
# S3 Bucket for upload
S3_BUCKET=''
SRA=''
UP_REGEX="*.[123].fq*"

# Script Arguments -DPU
UL_ARGS=''

# Output options -do
WORKDIR="$PWD"
OUTNAME="$SRA"

while getopts k:s:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    k)
      S3_BUCKET=$OPTARG
      ;;
    s)
      SRA=$OPTARG
      ;;
    U)
      UL_ARGS=$OPTARG
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

if [[ ( -z "$S3_BUCKET" ) ]]
then
  echo "-k S3_BUCKET is required."
  usage
fi

if [[ ( -z "$SRA" ) ]]
then
  echo "-s SRA/Accession id is required."
  usage
fi

# Copy all files matching UP_REGEX
# in current dir to s3
aws s3 cp --recursive \
  --exclude "*" \
  --include "$UP_REGEX" \
  ./  "$S3_BUCKET/$SRA"

