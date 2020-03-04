#!/bin/bash
# run_bowtie2
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
  echo "Usage: run_bowtie2.sh -1 READ1.fq(.gz) -2 READ2.fq -x GENOME [-LISPF] -o <output_prefix>"
  echo "   or: run_bowtie2.sh -0 READ0.fq(.gz) -x GENOME [-LISPF] -o <output_prefix>"
  echo ""
  echo "    Default behaviour is to retain mapped-reads and their unmapped-pairs only"
  echo ""
  echo "    Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -0 <0.fq>"
  echo "    -1    fastq paired-end reads 1"
  echo "    -2    fastq paired-end reads 2"
  echo "    -0    fastq unpaired reads"
  echo ""
  echo "    bowtie2 Alignment Parameters"
  echo "    -x    Name of bt2-indexed genome in <$WORKDIR>"
  echo "    -a    Aligner arguments  [--very-sensitive-local]"
  echo "    -p    N parallel threads [1]"
  echo ""
  echo "    Read Group (meta-data) required for GATK"
  echo "    -L    Library Name, Human readable (output_prefix)"
  echo "    -I    Group ID, SRA Accession Number"
  echo "    -S    Sample Name / Sample ID"
  echo "    -P    Population ID, or Study ID"
  echo "    -F    Platform ID [ILLUMINA]"
  echo ""
  echo "    Output options"
  echo "    -d    Working directory [pwd]"
  echo "    -o    <output_filename_prefix>"
  echo ""
  echo "    Outputs: <output_prefix>.bam"
  echo "             <output_prefix>.bam.bai"
  echo "             <output_prefix>.flagstat"
  echo ""
  echo "ex: ./run_bowtie2.sh -0 unpaired.fq -x hg38 -o testLib -I SRAX -S example -P silico"
  echo "ex: ./run_bowtie2.sh -1 toy.1.fq -2 toy.2.fq -x hgr1 -o toyLib -I SRAX -S example -P silico"
  exit 1
}

# PARSE INPUT =============================================
# Initialize input options -123
# Input Fastq files - paired 1+2 | unpaired 0
FQ1=""
FQ2=""
FQ0=""

# bowtie2 run parameters -xap
GENOME=""
BT2_ARG="--very-sensitive-local"
THREADS="1"

# Read Group meta-data (required for GATK) -LISPF
RGLB=""
RGID=""
RGSM=""
RGPO=""
RGPL="ILLUMINA"

# Output options -do
WORKDIR=''
OUTNAME=''

while getopts h0:1:2:x:a:p:L:I:S:P:F:d:o: FLAG; do
  case $FLAG in
    # Fastq Options ---------
    0)
      FQ0=$OPTARG
      ;;
    1)
      FQ1=$OPTARG
      ;;
    2)
      FQ2=$OPTARG
      ;;
    # bowtie2 options -------
    x)
      GENOME=$OPTARG
      ;;
    a)
      BT2_ARG=$OPTARG
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
# Logging to STDOUT
echo " -- bowtie2 Alignment Pipeline -- "
echo " date:    $(date)"
echo " version: $PIPE_VERSION "
echo " ami:     $AMI_VERSION  "
echo " output:  $OUTNAME"

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

echo " genome:  $GENOME"
echo " bt2 arguments:   $BT2_ARG -p $THREADS"
echo " Read Group LB:   $RGLB"
echo "            ID:   $RGID"
echo "            SM:   $RGSM"
echo "            PO:   $RGPO"
echo "            PL:   $RGPL"

echo""
echo 'Initializing ...'
echo ""
cd $WORKDIR

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

  echo "Alignment complete."

  # Clean-up FQ files to save space
  #rm $FQ1 $FQ2 

  # Extract aligned reads header
  samtools view -H aligned_unsorted.bam > align.header.tmp

  # Extract Mapped Reads and their unmapped pairs
    samtools view -b -F 4 aligned_unsorted.bam > align.F4.bam  # mapped
    samtools view -b -f 4 -F 8 aligned_unsorted.bam > align.f4F8.bam # unmapped-pair

  # Re-compile bam output
    samtools cat -h align.header.tmp -o align.tmp.bam align.F4.bam align.f4F8.bam
    samtools sort -@ $THREADS -O BAM align.tmp.bam > "$OUTNAME".bam

  # Flagstat and index
  samtools flagstat "$OUTNAME".bam > "$OUTNAME".flagstat
  samtools index "$OUTNAME".bam

  # OUTPUT: $OUTNAME.bam
  # OUTPUT: $OUTNAME.bam.bai
  # OUTPUT: $OUTNAME.flagstat


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

  echo "Alignment complete."

  # Clean-up FQ files to save space
  #rm $FQ0

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

fi

