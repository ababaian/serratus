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
  echo "Usage: run_merge.sh -b <bam-block prefix regex> -M '<samtools merge args>' [OPTIONS]"
  echo ""
  echo "    -h    Show this help/usage message"
  echo ""
  echo "    Required Parameters"
  echo "    -b    A regex string which matches the prefix the bam-block files in"
  echo"           working in base directory"
  echo "          i.e."
  echo "          if merging SRR1337.001.bam, SRR1337.002.bam, ... , SRR1337.666.bam"
  echo "          then  [-b 'SRR1337'] or  [-b 'SRR*'] will both work"
  echo "    -s    SRA Accession"
  echo ""
  echo "    Merge Parameters"
  echo "    -n    parallel CPU threads to use where applicable  [1]"
  echo "    -i    Flag. Do not generate bam.bai index file"
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

# Scheduler / Container Parameters
SCHED='localhost:8000'
S3_BUCKET='s3://serratus-public/out'

# Merge Options
THREADS='1'
INDEX='true'
FLAGSTAT='true'

# Inputs (manual overwrite)
SRA=''
S3_BAMS=''
GENOME=''

# AWS Options -ak
AWS_CONFIG='TRUE'

# Script Arguments -AU
MERGE_ARGS=''

# Output options -do
BASEDIR="/home/serratus"
OUTNAME="$SRA"

while getopts u:k:n:s:b:g:M:d:o:whif FLAG; do
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
      INDEX="false"
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

# Fastq required inputs: FQ1 and FQ2 or FQ0
paired_run=''
if [ -z "$FQ1" ] & [ -z "$FQ2" ]
then
  if [ -z "$FQ0" ]
  then
    echo "Fastq input required. Either -1 & -2, or -0"
    usage
  else
    paired_run="F"
  fi
else
  paired_run="T"
fi

# Required parameters
if [ -z "$GENOME" ]; then
  echo "-x <genome> required."
  usage
fi

if [[ -z "$RGID" || -z "$RGSM" || -z "$RGPO" ]]
then
  echo "-I -S -P read-groups required"
  usage
fi

if [ -z "$OUTNAME" ]
then
  echo "-o <output_name> required"
  usage
fi

# If RGLB is not defined; use <output_name>
if [ -z "$RGLB" ]
then
  RGLB=$OUTNAME
fi

# ALIGN ===================================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Logging to STDOUT
echo " -- bowtie2 Alignment Pipeline -- "
echo " date:    $(date)"
echo " version: $PIPE_VERSION "
echo " ami:     $AMI_VERSION  "
echo " output:  $OUTNAME"
echo " run-id:  $RUNID"

if [ $paired_run = "T" ]
then
  echo " type:    paired"
  echo " fq1:     $FQ1"
  echo " fq2:     $FQ2"
  echo " output:  $OUTNAME.bam"
else
  echo " type:    single-end"
  echo " fq0:     $FQ0"
  echo " output:  $OUTNAME.se.bam"
fi

echo""
echo 'Initializing ...'
echo ""

if [ $paired_run = "T" ]
then
  # Paired Read Alignment
  echo "bowtie2 $BT2_ARG -p $THREADS \\"
  echo "  --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \\"
  echo "  --rg PL:$RGPL --rg PU:$RGPU \\"
  echo "  -x hgr1 -1 $FQ1 -2 $FQ2 | \\"
  echo "samtools view -bS - > aligned_unsorted.bam"

  bowtie2 $BT2_ARG -p $THREADS \
    --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
    --rg PL:$RGPL --rg PU:$RGPU \
    -x hgr1 -1 $FQ1 -2 $FQ2 | \
    samtools view -bS - > aligned_unsorted.bam

  echo ""
  echo "Alignment complete."

  if [ $DEBUG = "F" ]
  then
    echo "clearing fq-files"
    # Clean-up FQ files to save space
    rm $FQ1 $FQ2
  fi

  echo "Extracting mapped reads + unmapped pairs"

  # Extract aligned reads header
    samtools view -H aligned_unsorted.bam > align.header.tmp

  # Extract Mapped Reads and their unmapped pairs
    samtools view -b -F 4 aligned_unsorted.bam > align.F4.bam  # mapped
    samtools view -b -f 4 -F 8 aligned_unsorted.bam > align.f4F8.bam # unmapped-pair

  # Re-compile bam output
    samtools cat -h align.header.tmp -o align.tmp.bam align.F4.bam align.f4F8.bam
    samtools sort -@ $THREADS -O BAM align.tmp.bam > "$OUTNAME".bam

  echo "Flagstat and index of output"

  # Flagstat and index
    samtools flagstat "$OUTNAME".bam > "$OUTNAME".flagstat
    samtools index "$OUTNAME".bam

  # OUTPUT: $OUTNAME.bam
  # OUTPUT: $OUTNAME.bam.bai
  # OUTPUT: $OUTNAME.flagstat

  if [ $DEBUG = "F" ]
  then
    # Clean-up temporary files
    rm *tmp aligned_unsorted.bam align.F4.bam align.f4F8.bam
  fi

else
  # Unpaired Read Alignment
  echo "bowtie2 $BT2_ARG -p $THREADS \\"
  echo "  --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \\"
  echo "  --rg PL:$RGPL --rg PU:$RGPU \\"
  echo "  -x hgr1 -0 $FQ0 | \\"
  echo "samtools view -bS - > aligned_unsorted.bam"

  bowtie2 $BT2_ARG -p $THREADS \
    --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
    --rg PL:$RGPL --rg PU:$RGPU \
    -x hgr1 -0 $FQ0 | \
    samtools view -bS - > aligned_unsorted.bam

  echo ""
  echo "Alignment complete."

  if [ $DEBUG = "F" ]
  then
    echo "clearing fq-files"
    # Clean-up FQ
    rm $FQ0
  fi
  
  echo "Extracting mapped reads + unmapped pairs"

  ## NOTE: Possibly skip sort / index steps if
  ## using fq/bam-blocks and not complete bam files
  # samtools flagstat > "$OUTNAME".flagstat
  # samtools view -bh -F 4 aligned_unsorted > "$OUTNAME".se.bam
  
  # Extract Mapped Reads and sort (flag 0x4)
  samtools view -bh -F 4 aligned_unsorted.bam | \
  samtools sort -@ $THREADS -O BAM - > "$OUTNAME".se.bam

  # Flagstat and index
  samtools flagstat aligned_unsorted.bam > "$OUTNAME".se.flagstat
  samtools index "$OUTNAME".se.bam

  # OUTPUT: $OUTNAME.se.bam
  # OUTPUT: $OUTNAME.se.bam.bai
  # OUTPUT: $OUTNAME.se.flagstat

  if [ $DEBUG = "F" ]
  then
    # Clean-up temporary files
    rm *tmp aligned_unsorted.bam
  fi
fi
