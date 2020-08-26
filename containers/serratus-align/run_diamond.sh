#!/bin/bash
# run_diamond
#
# Base: serratus-Aligner (0.1)
# AMI : ami-059b454759561d9f4
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.1.4"

#Usage function
function usage {
  echo ""
  echo "Usage: run_diamond.sh -1 READ_1.fq(.gz) -2 READ_2.fq -x <genome> -o <output_prefix>"
  echo "   or: run_bowtie2.sh -3 READS.fq(.gz) -x <genome>"
  echo ""
  echo "    Default behaviour is to retain mapped-reads only"
  echo "      paired-end data is not supported"
  echo ""
  echo "    Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -3 <Un.fq>"
  echo "    -1    path to fastq paired-end reads 1"
  echo "    -2    path to fastq paired-end reads 2"
  echo "    -3    path to fastq unpaired reads"
  echo ""
  echo "    bowtie2 Alignment Parameters"
  echo "    -x    path to _name_ of dmnd-indexed genome (exclude .dmnd extension)"
  echo "    -a    Aligner arguments  [not implemented]"
  echo "    -p    N parallel threads [1]"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo "    -!    Debug mode, will not rm intermediate files"
  echo ""
  echo "    Outputs: <output_prefix>.bam"
  echo ""
  exit 1
}

# PARSE INPUT =============================================
# Initialize input options -123
# Input Fastq files - paired 1+2 | unpaired 0
FQ1=""
FQ2=""
FQ3=""

# bowtie2 run parameters -xap
GENOME=""
THREADS="1"

# Read Group meta-data (required for GATK) -LISPF
RGLB=""
RGID=""
RGSM=""
RGPO=""
RGPL="ILLUMINA"

# Output options -do
WORKDIR="$PWD"
OUTNAME=''
DEBUG='F'

while getopts h3:1:2:x:a:p:L:I:S:P:F:d:o:! FLAG; do
  case $FLAG in
    # Fastq Options ---------
    3)
      FQ3=$(readlink -f $OPTARG)
      ;;
    1)
      FQ1=$(readlink -f $OPTARG)
      ;;
    2)
      FQ2=$(readlink -f $OPTARG)
      ;;
    # diamond options -------
    x)
      GENOME=$OPTARG
      ;;
    a)
      if [ ! -z "$OPTARG" ]
      then
        # Use user-input 
        ALIGN_ARG=$OPTARG
      fi
      ;;
    p)
      THREADS=$OPTARG
      ;;
    # Read Groups (not used)
    L)
      RGLB=$OPTARG
      ;;
    I)
      RGID=$OPTARG
      ;;
    S)
      RGSM=$OPTARG
      ;;
    P)
      RGPO=$OPTARG
      ;;
    F)
      RGPL=$OPTARG
      ;;
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

# Fastq required inputs: FQ1 and FQ2 or FQ3
paired_run=''
if [ -z "$FQ1" ] & [ -z "$FQ2" ]
then
  if [ -z "$FQ3" ]
  then
    echo "Fastq input required. Either -1 & -2, or -U"
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

if [ -z "$OUTNAME" ]
then
  echo "-o <output_name> required"
  usage
fi

# ALIGN ===================================================
# Unit test:
# yum install -y tar gzip
# aws s3 cp s3://serratus-public/test-data/fq/ERR2756788.fq ./
#
# aws s3 cp s3://serratus-public/rce/mprot/tarball/mprot.tz ./
# cp out/mall.fa ./
## Make diamond database (dmnb)
# diamond makedb --in mall.fa -d mall
#
# mkdir tmp
# cat ERR2756788.fq |\
# diamond blastx \
#   -d mall.dmnd \
#   --unal 0 \
#   -k 1 \
#   -p 1 \
#   -b 0.2 \
#   -f 6 qseqid sseqid qstart qend qlen sstart send slen pident evalue btop cigar qstrand qseq sseq \
#   > tmp.bam
## 3357 Alignments
## 7cd7eefa94d29ef9f1e3d588dc79713b  tmp.bam
# -----------------

# Data is labeled as "bam" but it is not
# data is saved as a diamond TSV (BLAST)
# -k max hits
# -p threads
# -b block-size

cat *.fq* |\
diamond blastx \
  -d "$GENOME".dmnd \ 
  --unal 0 \
  -k 1 \
  -p 1 \
  -b 0.4 \
  -f 6 qseqid sseqid qstart qend qlen sstart send slen pident evalue btop cigar qstrand qseq sseq \
  > "$OUTNAME".bam

  # --unal 0 Do not report unmapped reads
  # -k 1 Maximum hits per query sequence. Reduce output file size for Cov+.
  # -p 1 Single-threaded for t2.micro or t2.nano .
  # -b 0.2 Limits memory, RAM used is roughly <6×b in Gb, ~1.2 Gb per core
  # -t /tmp Temp directory. Required to avoid bug with input from /dev/stdin.
  # -f 6 <field> BLAST tabular output format
