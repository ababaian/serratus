#!/bin/bash
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

#Usage function
function usage {
  echo ""
  echo "Usage: docker exec serratus-dl -s <SRA Accession> [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    SRA Accession for sratools"
  echo "    -s    SRA Accession for downloading SRA-archive file"
  echo ""
  echo "    AWS / S3 Bucket parameters"
  echo "    -a    Flag. Imports AWS IAM from this ec2-instance for container"
  echo "          EC2 instance must have been launched with correct IAM"
  echo "          (No alternative yet, hard set to TRUE)"
  echo "    -k    S3 bucket name to store fq-blocks [s3://serratus-public]"
  echo ""
  echo "    -D    String of arguments to pass to 'run_download.sh <ARG>'"
  echo "    -P    String of arguments to pass to 'run_split.sh <ARG>'"
  echo "    -U    String of arguments to pass to 'run_upload.sh <ARG>'"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory in container [/home/serratus/]"
  echo "    -o    <output_filename_prefix> [Defaults to SRA_ACCESSION]"
  echo ""
  echo "    Outputs Uploaded to s3: "
  echo "          <output_prefix>.fq.xxxx ... <output_prefix.fq.yyyy"
  echo ""
  echo "ex: docker exec serratus-dl -S SRR1337"
  exit 1
}

# PARSE INPUT =============================================
# SRA Accession -S
SRA=''

# AWS Options -ak
AWS_CONFIG='TRUE'
S3_BUCKET=''

# Script Arguments -DPU
DL_ARGS=''
SPLIT_ARGS=''
UL_ARGS=''

# Output options -do
WORKDIR="/home/serratus"
OUTNAME="$SRA"

while getopts s:ak:D:P:U:d:oh FLAG; do
  case $FLAG in
    # Input Options ---------
    s)
      SRA=$(readlink -f $OPTARG)
      ;;
    a)
      AWS_CONFIG='TRUE'
      ;;
    k)
      S3_BUCKET=$(readlink -f $OPTARG)
      ;;
    # SCRIPT ARGUMENTS -------
    D)
      DL_ARGS=$(readlink -f $OPTARG)
      ;;
    P)
      SPLIT_ARGS=$OPTARG
      ;;
    U)
      UL_ARGS='TRUE'
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

echo "============================"
echo "    serratus-dl Pipeline    "
echo "============================"
echo " date:      $(date)"
echo " version:   $PIPE_VERSION"
echo " ami:       $AMI_VERSION"
echo " container: $CONTAINER_VERSION"
echo " sra:       $SRA"
echo " args:      ..."
echo "$@"
echo ""

# Default /home/serratus
cd $WORKDIR

# AUTHENTICATE AWS ========================================
echo "  Authenticating AWS credentials"

# Create aws key (read from ec2 IAM) (run at start-up of container)
AWS_ACCESSKEYID=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.AccessKeyId' )
AWS_SECRETKEY=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.SecretAccessKey' )
echo User name,Password,Access key ID,Secret access key,Console login link > key.csv
echo default,,$AWS_ACCESSKEYID,$AWS_SECRETKEY, >> $WORKDIR/key.csv
chmod 400 $WORKDIR/key.csv

# Pass credentials to sra-toolkit via vdb-config
# Current version requires a manual
# vdb-config -i initialization
# ugly hack is to copy a blank config file and bypass this
mkdir -p /root/.ncbi
aws s3 cp s3://serratus-public/VDB_user-settings.mkfg /root/.ncbi/user-settings.mkfg
chmod 500 /root/.ncbi/user-settings.mkfg
vdb-config --accept-aws-charges yes \
  --report-cloud-identity yes \
  --set-aws-credentials $WORKDIR/key.csv

# Download AWS S3 test token
aws s3 cp s3://serratus-public/aws-test-token.jpg $WORKDIR/aws-test-token.jpg

if [ ! -s './aws-test-token.jpg' ]
then
  echo "    ERROR: AWS Test did not download"
  echo "    Ensure the EC2 instance has correct IAM permissions"
  echo "      - requires S3 Read/Write"
  usage
fi

echo '  ...authentication token was download successfully'
echo ""

# RUN DOWNLOAD ============================================
echo "  Running -- run_download.sh --"
echo "  ./scripts/run_download.sh -s $SRA $DL_ARGS"

./scripts/run_download.sh -s $SRA $DL_ARGS

# RUN SPLIT ===============================================
# Add FQ0 vs. FQ1+FQ2 logic here
echo "  Running -- run_split.sh --"
echo "  ./scripts/run_split.sh -o $OUTNAME $SPLIT_ARGS"

# RUN UPLOAD ==============================================
echo "  Running -- run_upload.sh --"
echo "  ./scripts/run_download.sh -k $S3_BUCKET $UL_ARGS"

# CLEAN-UP ================================================
#cd $WORKDIR; rm key.csv aws-test-token.jpg