#!/bin/bash
# run_merge
#
# Base: serratus-Aligner (>0.1)
# AMI : ami-059b454759561d9f4
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.1"
AMI_VERSION='ami-059b454759561d9f4'

function usage {
  echo ""
  echo "Usage: run_merge.sh -b <bam-block regex> -M '<samtools merge args>' [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Parameters"
  echo "    -s    SRA Accession / output name"
  echo ""
  echo "    Merge Parameters"
  echo "    -b    A regex string which matches the prefix the bam-block files in"
  echo"           working in base directory [ *.bam ]"
  echo "          i.e."
  echo "          if merging SRR1337.001.bam, SRR1337.002.bam, ... , SRR1337.666.bam"
  echo "          then  [-b 'SRR1337*'] or  [-b 'SRR*'] will both work"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo "    -i    Flag. Generate bam.bai index file"
  echo "    -f    Flag. Do not generate flagstat summary file"
  echo "    -r    Flag. Do not chromosome sort output"
  echo ""
  echo "    Arguments as string to pass to `samtools merge` command"
  echo "    -M    String of arguments"
  echo ""
  echo "    Output options"
  echo "    -d    Work directory in container with bam-block [pwd]"
  echo "    -o    <output_filename_prefix> [Defaults to SRA_ACCESSION]"
  echo ""
  echo "    Outputs a sorted Uploaded to s3: "
  echo "          <output_prefix>.bam, <output_prefix>.bam.bai, <output_prefix>.flagstat"
  echo ""
  echo "ex: bash run_merge.sh -b 'SRR*' -M '-l 9'"
  exit 1
}


# PARSE INPUT =============================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Run Parameters
BAMREGEX='*.bam'
SRA=''

# Merge Options
THREADS='1'
INDEX='false'
FLAGSTAT='true'
SORT='true'

# Script Arguments -M
MERGE_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts b:s:nifrM:d:o:h FLAG; do
  case $FLAG in
    # Scheduler Options -----
    u)
      SCHED=$OPTARG
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
      S3_BAMS=$OPTARG
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

# Check inputs --------------
# Required parameters
if [ -z "$SRA" ]
then
  echo "-s <SRA / Outname> required"
  usage
fi


# MERGE ===================================================
# Final Output Bam Files
OUTBAM="$SRA.bam"

# Create tmp file of bam files to merge
ls $BAMREGEX > bam.list

# Initial merge gives un-sorted; with ugly header
samtools merge -n -@ $THREADS \
  -b bam.list \
  temporary.bam

# Extract header from first file
samtools view -H $(head -n1 bam.list) > 0.header.sam

if [[ "$SORT" -eq 'true' ]]
then
  # Sort + Reheader
  samtools view -h temporary.bam | \
  python3 /home/serratus/serratus_summarizer.py /dev/stdin $SRA.summary /dev/stdout | \
  samtools sort -@ $THREADS - | \
  samtools reheader 0.header.sam - >\
  merged_sorted.bam

  mv merged_sorted.bam temporary.bam
  rm 0.header.sam
else
  # Re-header only
  samtools view -h temporary.bam | \
  python3 /home/serratus/serratus_summarizer.py /dev/stdin $SRA.summary /dev/stdout | \
  samtools reheader 0.header.sam - > \
  reheader.bam

  mv reheader.bam temporary.bam
  rm 0.header.sam
fi

# Set output name
mv temporary.bam $OUTBAM

# Clear bam-blocks
xargs rm -r <bam.list; rm bam.list

# Only final bam file in current directory

if [[ "$INDEX" -eq 'true' ]]
then
  samtools index $OUTBAM
fi

if [[ "$FLAGSTAT" -eq 'true' ]]
then
  samtools flagstat $OUTBAM > $SRA.flagstat

  echo "  Post Merge Flagstat"
  cat $SRA.flagstat
  echo ''
fi

# end of script
