#!/bin/bash
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
  echo "Usage: run_mblast-index.sh -x GENOME.fa -o <output_prefix>"
  echo ""
  echo "    Magic-blast alignment Parameters"
  echo "    -x    path to _name_ of genome (exclude .fa extension)"
  echo "    -i    Index arguments  []"
  echo "    -p    N parallel threads [1]"
  echo ""
  echo "    Output options"
  echo "    -r    Do not create REVERSE control index [false]"
  echo "          Reverse, non-compliment input sequence."
  echo "          Headers prefixed with 'REVERSE_'"
  echo "          Output appended with '.r' suffix"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo ""
  echo "    Outputs: "
  echo ""
  echo "ex: ./run_mblast-index.sh -x seq/cov0.fa -r -o cov0 "
  exit 1
}

# PARSE INPUT =============================================
# Initialize input options -123

# bowtie2 run parameters -xap
GENOME=""
IDX_ARG=""
THREADS="1"

# Output options -do
ADD_REVERSE='true'
WORKDIR="$PWD"
OUTNAME=''
DEBUG='F'

while getopts hx:i:p:r:d:o: FLAG; do
  case $FLAG in
    # bowtie2 options -------
    x)
      GENOME=$(readlink -f $OPTARG)
      ;;
    i)
      IDX_ARG=$OPTARG
      ;;
    p)
      THREADS=$OPTARG
      ;;
    # Output opt. -----------
    r)
      # do not create if flag set
      REVERSE='false'
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
# Required parameters

if [ -z "$GENOME" ]
then
  echo "Genome fasta input required. -x"
  usage
fi

if [ -z "$OUTNAME" ]
then
  echo "-o <output_name> required"
  usage
fi

if [ "$ADD_REVERSE" = "true" ]
then
  TMP="$OUTPUT".r
  OUTPUT=$TMP
else
  TMP="$OUTPUT"
  OUTPUT=$TMP
fi

# ALIGN ===================================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Logging to STDOUT
echo " -- mblast-index Pipeline -- "
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
  echo " output:  $OUTNAME.bam"
fi

# ---------------------------

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
  echo "  -x hgr1 -U $FQ0 | \\"
  echo "samtools view -bS - > aligned_unsorted.bam"

  bowtie2 $BT2_ARG -p $THREADS \
    --rg-id $RGID --rg LB:$RGLB --rg SM:$RGSM \
    --rg PL:$RGPL --rg PU:$RGPU \
    -x hgr1 -U $FQ0 | \
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
  # samtools view -bh -F 4 aligned_unsorted > "$OUTNAME".bam
  
  # Extract Mapped Reads and sort (flag 0x4)
  samtools view -bh -F 4 aligned_unsorted.bam | \
  samtools sort -@ $THREADS -O BAM - > "$OUTNAME".bam

  # Flagstat and index
  samtools flagstat aligned_unsorted.bam > "$OUTNAME".flagstat
  samtools index "$OUTNAME".bam

  # OUTPUT: $OUTNAME.bam
  # OUTPUT: $OUTNAME.bam.bai
  # OUTPUT: $OUTNAME.flagstat

  if [ $DEBUG = "F" ]
  then
    # Clean-up temporary files
    rm aligned_unsorted.bam
  fi
fi
