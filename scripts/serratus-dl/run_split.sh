#!/bin/bash
# run_split
#
set -eu
# Base: serratus-Downloader (>0.1.1)
# Amazon Linux 2 with Docker
# AMI : aami-0fdf24f2ce3c33243
# Container: 
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.1"
AMI_VERSION='ami-0fdf24f2ce3c33243'

#Usage function
function usage {
  echo ""
  echo "Usage: run_split.sh -b <input.bam> -o <output_prefix> [OPTIONS]"
  echo "   or: run_split.sh -f <unpaired.fq>  -o <output_prefix> [OPTIONS]"
  echo "   or: run_split.sh -1 <input.1.fq> -2 <input.2.fq>  -o <output_prefix> [OPTIONS]"
  echo ""
  echo "    BAM input req"
  echo "    -b    path to input bam file. Auto-detect single/paired-end"
  echo ""
  echo "    Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -f <in.fq>"
  echo "    -f    path to single-end reads"
  echo "    -1    path to fastq paired-end reads 1"
  echo "    -2    path to fastq paired-end reads 2"
  echo ""
  echo "    fq-block parameters"
  echo "    -n    reads per fq-block [1000000]"
  echo "          approx 2.2Mb per 10k reads (single-end)"
  echo "              or 220Mb per 1M  reads"
  echo "    -p    N parallel CPU threads [1]"
  echo "    -z    flag to NOT use gzip fq-blocks [ Default is .gz ]"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo "    (N/A)-!    Debug mode, will not rm intermediate files"
  echo ""
  echo "    Outputs: <output_prefix>.bam"
  echo "             <output_prefix>.bam.bai"
  echo "             <output_prefix>.flagstat"
  echo ""
  echo "ex: ./run_split.sh -b tester.bam -o testFromBam"
  echo "ex: ./run_split.sh -f SRA1337.fq -n 10000 -p 8 -z -o SRA1337"
  exit 1
}

# PARSE INPUT =============================================
# Initialize input options -b123
BAM=""
# Input Fastq files - paired 1+2 | unknown 0
FQ1=""
FQ2=""
FQ3=""

# fq-block parameters -nzp
BLOCKSIZE=1000000
GZIP_FLAG="true"
THREADS="1"

# Output options -do
WORKDIR="$PWD"
OUTNAME=''
DEBUG='F'

while getopts b:f:1:2:n:p:zd:o:!hz FLAG; do
  case $FLAG in
    # Input Options ---------
    b)
      BAM=$(readlink -f $OPTARG)
      ;;
    f)
      FQ3=$(readlink -f $OPTARG)
      ;;
    1)
      FQ1=$(readlink -f $OPTARG)
      ;;
    2)
      FQ2=$(readlink -f $OPTARG)
      ;;
    # fq-block options -------
    n)
      BLOCKSIZE=$OPTARG
      ;;
    z)
      GZIP_FLAG='false'
      ;;
    p)
      THREADS=$OPTARG
      ;;
    # output options -------
    d)
      WORKDIR=$OPTARG
      ;;
    o)
      OUTNAME=$OPTARG
      ;;
    !)
      DEBUG="T"
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

if [[ ( -z "$BAM" ) && ( -z "$FQ3" ) && ( -z "$FQ1" || -z "$FQ2" ) ]]
then
  echo "(-b) .bam         OR"
  echo "(-f) .fq unpaired OR"
  echo "(-1 -2) 1.fq and 2.fq paired-end fastq input required"
  usage
fi

if [ -z $OUTNAME ]
then
  echo "(-o) output prefix is required."
  usage
fi

# SPLIT ===================================================

# Logging to STDOUT
echo " -- fq-split Alignment Pipeline -- "
echo " date:      $(date)"
echo " version:   $PIPE_VERSION"
echo " ami:       $AMI_VERSION"
echo " workdir:   $WORKDIR"
echo " input(s):  $BAM $FQ3 $FQ1 $FQ2"
echo " output:    $OUTNAME"
echo " blockSize: $BLOCKSIZE"

echo""
echo 'Initializing ...'
echo ""

# fq-block generation function ----------------------------
  fq-block-generate () {
    # Input a fastq file; output fq-blocks
    echo "  Spliting $1 into fq-blocks of $BLOCKSIZE reads."

    # Split in-fq to $BLOCKSIZE reads per file (4lines/read)
    # will generate n * input.fq.abcdefghi.gz files
    let LINESIZE=4*BLOCKSIZE
    split -d -a 10 -l $LINESIZE $1 "$1".
    rm $1

    # gzip fq-blocks in parallel
    if [ $GZIP_FLAG = "true" ]
    then
      pigz -p $THREADS "$1"*
    fi
  }

# BAM-based input -----------------------------------------
if [[ -s $BAM && -n $BAM ]]
then
  echo "  Processing bam-input"
  # Set created FQ filename variables
  FQ1=$OUTNAME.1.fq
  FQ2=$OUTNAME.2.fq
  FQ3=$OUTNAME.3.fq

  # Sort input bam by read-name (for paired-end)
  samtools sort -n -l 0 $BAM |
  samtools fastq \
    -0 $FQ3 \
    -1 $FQ1 -2 $FQ2 \
    -@ $THREADS -

  # Count how many lines in each file
  lc1=$(wc -l $FQ1 | cut -f1 -d' ' -)
  lc2=$(wc -l $FQ2 | cut -f1 -d' ' -)
  lc3=$(wc -l $FQ3 | cut -f1 -d' ' -)

  # For Paired-End Reads; ensure equal read-pairs
  # to continue
  if [ $lc1 != $lc2 ]
  then
    echo "Early Error: number of lines in 1.fq != 2.fq"
    echo "There are unqueal paired-end reads or mixed paired/single"
    echo "Force input as unpaired or convert to ubam and retry"
    exit 2
  else
    echo "all good!"
  fi

  if [ $lc3 != "0" ]
  then
    fq-block-generate $FQ3
  else
    echo "$FQ3 is empty... skipping"
  fi

  if [ $lc1 != "0" ]
  then
    fq-block-generate $FQ1
  else
    echo "$FQ1 is empty... skipping"
  fi

  if [ $lc2 != "0" ]
  then
    fq-block-generate $FQ2
  else
    echo "$FQ2 is empty... skipping"
  fi

# Paired-End FQ input -------------------------------------
elif [[ ( -s $FQ1 && -n $FQ1 ) && ( -s $FQ2 && -n $FQ2 ) ]]
then
  echo "  Processing paired-end fq-input"
  # Link fq to RUNID folder
  ln -s $FQ1 $OUTNAME.1.fq
  ln -s $FQ2 $OUTNAME.2.fq
  FQ1=$OUTNAME.1.fq
  FQ2=$OUTNAME.2.fq

  # Count how many lines in each file
  lc1=$(wc -l $FQ1 | cut -f1 -d' ' -)
  lc2=$(wc -l $FQ2 | cut -f1 -d' ' -)

  # For Paired-End Reads; ensure equal read-pairs
  # to continue
  if [ $lc1 != $lc2 ]
  then
    echo "Early Error: number of lines in 1.fq != 2.fq"
    echo "There is a loss of paired-end reads."
    echo "Force-unpaired or recompile"
    exit 2
  else
    echo "all good!"
  fi

  if [ $lc1 != "0" ]
  then
    fq-block-generate $FQ1
  else
    echo "$FQ1 is empty... skipping"
  fi

  if [ $lc2 != "0" ]
  then
    fq-block-generate $FQ2
  else
    echo "$FQ2 is empty... skipping"
  fi

# Single-End FQ input -------------------------------------
elif [[ -s $FQ3 && -n $FQ3 ]]
then
  # Link fq to RUNID folder
  ln -s $FQ3 $OUTNAME.3.fq
  FQ3=$OUTNAME.3.fq

  echo "  Processing single-end fq-input"
  # Count how many lines in each file
  lc3=$(wc -l $FQ3 | cut -f1 -d' ' -)

  # Assume a single fastq file is single-end reads
  if [ $lc3 != "0" ]
  then
    fq-block-generate $FQ3
  else
    echo "$FQ3 is empty... skipping"
  fi
else
  echo "Error: input(s) do not exist or are non-readable"
  # Something is wrong, one of BAM, FQ3, FQ1 and FQ2 should exist.
  exit 3
fi

# :)

# Read paired-end FASTQ files into a BAM-stream, n-sort for splitting

#gatk FastqToSam -F1 NA12878.1.fq.gz -F2 NA12878.2.fq.gz -O test.bam -SM test_sample
