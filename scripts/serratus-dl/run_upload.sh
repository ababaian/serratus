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
  echo "    -k    S3 URL to upload data to [s3://serratus-public/fq-block]"
  echo "    -s    SRA Accession name"
  echo ""
  echo "    Arguments from serratus-dl"
  echo "    <-U>   Arguments as string from running serratus-dl passed to"
  echo "           this script (optional)"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [/home/serratus]"
  echo "    -o    <output_prefix>"
  echo ""
  echo "ex: ./run_upload.sh -k s3://serratus-public/fq-block -s SRA1337"
  echo ""
  echo "    Uploads to: s3:"
  echo ""
  exit 1
}

# PARSE INPUT =============================================
# S3 Bucket for upload
S3_BUCKET=''

# Script Arguments -DPU
UL_ARGS=''

# Output options -do
WORKDIR="/home/serratus"
OUTNAME="$SRA"

while getopts k:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    k)
      S3_BUCKET=$(readlink -f $OPTARG)
      ;;
    U)
      UL_ARGS=$(readlink -f $OPTARG)
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

echo " -- fq-split Alignment Pipeline -- "
echo " date:      $(date)"
echo " version:   $PIPE_VERSION"
echo " ami:       $AMI_VERSION"
echo " container: $CONTAINER_VERSION"
echo " sra:       $S3_BUCKET"
echo " args:      ..."
echo "$@"
echo ""

# Default /home/serratus
cd $WORKDIR
