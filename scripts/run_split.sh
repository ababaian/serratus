#!/bin/bash
# run_split
#
# Base: serratus-Aligner (>0.1)
# AMI : ami-059b454759561d9f4
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.0"
AMI_VERSION='ami-059b454759561d9f4'

#Usage function
function usage {
  echo ""
  echo "Usage: run_split.sh -b <input.bam> -o <output_prefix> [OPTIONS]"
  echo "   or: (N/A) run_split.sh -f <input.fq>  -n <int-reads-per-block>"
  echo "   or: (N/A) run_split.sh -1 <input.1.fq> -2 <input.2.fq>  -n <int-reads-per-block>"
  echo ""
  echo "    BAM input req"
  echo "    -b    path to input bam file. Auto-detect single/paired-end"
  echo ""
  echo "    (N/A) Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -f <in.fq>"
  echo "    -f    path to single-end or interleaved fastq file"
  echo "    -1    path to fastq paired-end reads 1"
  echo "    -2    path to fastq paired-end reads 2"
  echo ""
  echo "    fq-block parameters"
  echo "    -n    reads per fq-block [1000000]"
  echo "          approx 2.2Mb per 10k reads (single-end)"
  echo "              or 220Mb per 1M  reads"
  echo "    -p    N parallel threads [1]"
  echo "    -z    flag to gzip fq-blocks [F]"
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
  echo "ex: ./run_split.sh -b tester.bam -o testReads"
  echo "ex: ./run_split.sh -b testing.bam -n 10000 -p 8 -o tmp"
  exit 1
}

# PARSE INPUT =============================================
# Initialize input options -b012
BAM=""
# Input Fastq files - paired 1+2 | unknown 0
FQ1=""
FQ2=""
FQ0=""

# fq-block parameters -nzp
BLOCKSIZE=1000000
GZIP_FLAG="FALSE"
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
      FQ0=$(readlink -f $OPTARG)
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
      GZIP_FLAG='TRUE'
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

if [ -z $BAM ]
then
  echo "(-b) .bam input required for initial version."
  usage
fi

if [ -z $OUTNAME ]
then
  echo "(-o) output prefix is required."
  usage
fi

# SPLIT ===================================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Logging to STDOUT
echo " -- fq-split Alignment Pipeline -- "
echo " date:      $(date)"
echo " version:   $PIPE_VERSION"
echo " ami:       $AMI_VERSION"
echo " run-id:    $RUNID"
echo " input:     $BAM"
echo " output:    $OUTNAME"
echo " blockSize: $BLOCKSIZE"

echo""
echo 'Initializing ...'
echo ""


# Create RUNID folder in WORKDIR
mkdir -p $WORKDIR/$RUNID
cd $WORKDIR/$RUNID
ln -s $BAM ./

  # Sort input bam by read-name (for paired-end)
  samtools sort -n -l 0 $BAM |
  samtools fastq \
    -0 $OUTNAME.0.fq \
    -1 $OUTNAME.1.fq -2 $OUTNAME.2.fq \
    -@ $THREADS -

  # Count how many lines in each file
  lc0=$(wc -l $OUTNAME.0.fq | cut -f1 -d' ' -)
  lc1=$(wc -l $OUTNAME.1.fq | cut -f1 -d' ' -)
  lc2=$(wc -l $OUTNAME.2.fq | cut -f1 -d' ' -)

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

  fq-block-generate () {
      echo " Spliting $1 into fq-blocks of $BLOCKSIZE reads."

      # Split in-fq to $BLOCKSIZE reads per file (4lines/read)
      # will generate n * input.fq.abcdefghi.gz files
      let LINESIZE=4*BLOCKSIZE
      split -a 10 -l $LINESIZE $1 "$1".
      rm $1

      # gzip fq-blocks in parallel
      if [ $GZIP_FLAG = "TRUE" ]
      then
        pigz -n $THREADS "$1"*
      fi
  }

  if [ $lc0 != "0" ]
  then
    fq-block-generate $OUTNAME.0.fq
  else
    echo "$OUTNAME.0.fq is empty... skipping"
  fi

  if [ $lc1 != "0" ]
  then
    fq-block-generate $OUTNAME.1.fq
  else
    echo "$OUTNAME.1.fq is empty... skipping"
  fi

  if [ $lc2 != "0" ]
  then
    fq-block-generate $OUTNAME.2.fq
  else
    echo "$OUTNAME.2.fq is empty... skipping"
  fi


# :)