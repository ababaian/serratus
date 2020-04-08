#!/bin/bash
set -eux
# run_mblast
#
# Base: serratus-Aligner (>0.2)
# AMI : ami-059b454759561d9f4
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.2"
AMI_VERSION='ami-059b454759561d9f4'

#Usage function
function usage {
  echo ""
  echo "Usage: run_mblast.sh -1 READ1.fq(.gz) -2 READ2.fq -x GENOME -o <output_prefix>"
  echo "   or: run_mblast.sh -3 READ0.fq(.gz) -x GENOME -o <output_prefix>"
  echo ""
  echo "    Default behaviour is to retain mapped-reads and their unmapped-pairs only"
  echo ""
  echo "    Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -U <0.fq>"
  echo "    -1    path to fastq paired-end reads 1"
  echo "    -2    path to fastq paired-end reads 2"
  echo "    -3    path to fastq unpaired reads"
  echo ""
  echo "    magic-blast Alignment Parameters"
  echo "    -x    path to _name_ of genome (exclude .fa extension)"
  echo "    -a    Aligner arguments"
  echo "          [-splice F -no_unaligned -max_db_word_count 1000000]"
  echo "    -p    N parallel threads [1]"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo "    -!    Debug mode, will not rm intermediate files"
  echo ""
  echo "    Outputs: <output_prefix>.bam"
  echo "             <output_prefix>.bam.bai"
  echo "             <output_prefix>.flagstat"
  echo ""
  echo "ex: ./run_mblast.sh -3 ~/unpaired.fq -x ~/hg38 -o testLib -I SRAX -S example -P silico"
  echo "ex: ./run_mblast.sh -1 /scratch/toy.1.fq -2 /scratch/toy.2.fq -x /tmp/hgr1 -o toyLib -I SRAX -S example -P silico"
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
ALG_ARG="-splice F -no_unaligned -max_db_word_count 1000000"
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
    # bowtie2 options -------
    x)
      GENOME=$OPTARG
      ;;
    a)
      ALG_ARG=$OPTARG
      ;;
    p)
      THREADS=$OPTARG
      ;;
    # Read Groups -----------
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
echo " -- magic-blast Alignment Pipeline -- "
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
  echo " fq3:     $FQ3"
  echo " output:  $OUTNAME.bam"
fi

# ---------------------------

echo " genome:  $GENOME"
echo " mblast arguments: $ALG_ARG -p $THREADS"
echo " Read Group LB:   $RGLB"
echo "            ID:   $RGID"
echo "            SM:   $RGSM"
echo "            PO:   $RGPO"
echo "            PL:   $RGPL"

echo""
echo 'Initializing ...'
echo ""

if [ $paired_run = "T" ]
then
  ALG_ARG="-splice F -max_db_word_count 1000000"
  GENOME='cov0r'

  echo "magicblast -infmt fastq -paired $ALG_ARG "
  echo "  -query $FQ1 -query_mate $FQ2 "
  echo "  -db $GENOME | "
  echo "  samtools view -b - > aligned_unsorted.bam"

  magicblast -infmt fastq -paired $ALG_ARG \
    -query $FQ1 -query_mate $FQ2 \
    -db $GENOME | \
    samtools view -b -G 12 - > "$OUTNAME".bam

  echo ""
  echo "Alignment complete."

  if [ $DEBUG = "F" ]
  then
    echo "clearing fq-files"
    # Clean-up FQ files to save space
    rm $FQ1 $FQ2
  fi

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
    rm -f *tmp aligned_unsorted.bam align.F4.bam align.f4F8.bam
  fi

else
  # Paired is not true
  ALG_ARG="-splice F -max_db_word_count 1000000"
  GENOME='cov0r'

  echo "magicblast -infmt fastq $ALG_ARG "
  echo "  -query $FQ3 "
  echo "  -db $GENOME | "
  echo "  samtools view -b - > aligned_unsorted.bam"

  magicblast -infmt fastq -no_unaligned $ALG_ARG \
    -query $FQ3 \
    -db $GENOME | \
    samtools view -b - > "$OUTNAME".bam

  echo ""
  echo "Alignment complete."

    if [ $DEBUG = "F" ]
    then
      echo "clearing fq-files"
      # Clean-up FQ
      rm -f $FQ3
    fi
  
  echo "Extracting mapped reads + unmapped pairs"

  ## NOTE: Possibly skip sort / index steps if
  ## using fq/bam-blocks and not complete bam files
  # samtools flagstat > "$OUTNAME".flagstat
  # samtools view -bh -F 4 aligned_unsorted > "$OUTNAME".bam
  
  # Extract Mapped Reads and sort (flag 0x4)
  # samtools view -bh -F 4 aligned_unsorted.bam | \
  #samtools sort -@ $THREADS -O BAM - > "$OUTNAME".bam

  # Flagstat and index
  samtools flagstat "$OUTNAME
  ".bam > "$OUTNAME".flagstat
  samtools index "$OUTNAME".bam

  # OUTPUT: $OUTNAME.bam
  # OUTPUT: $OUTNAME.bam.bai
  # OUTPUT: $OUTNAME.flagstat

  if [ $DEBUG = "F" ]
  then
    # Clean-up temporary files
    rm -f aligned_unsorted.bam
  fi
fi
