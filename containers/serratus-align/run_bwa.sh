#!/bin/bash
# run_bowtie2
#
# Base: serratus-Aligner (>0.1)
# AMI : ami-059b454759561d9f4
# login: ec2-user@<ipv4>
# base: 9 Gb
#
PIPE_VERSION="0.1.3"

#Usage function
function usage {
  echo ""
  echo "Usage: run_bwa.sh -1 READ_1.fq(.gz) -2 READ_2.fq -x <genome> -L <RGLB> -I <RGID> -S <RGSM> -P <RGPO> -o <output_prefix>"
  echo "   or: run_bwa.sh -3 READS.fq(.gz) -x <genome> [-LISPF] -o <output_prefix>"
  echo ""
  echo "    Default behaviour is to retain mapped-reads and their unmapped-pairs only"
  echo ""
  echo "    Fastq input req: (-1 <1.fq> -2 <2.fq> ) || -3 <Un.fq>"
  echo "    -1    path to fastq paired-end reads 1"
  echo "    -2    path to fastq paired-end reads 2"
  echo "    -3    path to fastq unpaired reads"
  echo ""
  echo "    bwa-mem Alignment Parameters"
  echo "    -x    path to _name_ of bwa-indexed genome (exclude .fa extension)"
  echo "    -a    Aligner arguments  [-k 14]"
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
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo "    -!    Debug mode, will not rm intermediate files"
  echo ""
  echo "    Outputs: <output_prefix>.bam"
  echo ""
  echo "ex: ./run_bwa.sh -3 ~/unpaired.fq -x cov1r -o testLib -I SRAX -S example -P silico"
  echo "ex: ./run_bwa.sh -1 /scratch/toy.1.fq -2 /scratch/toy.2.fq -x ~/tmp/hgr1 -o toyLib -I SRAX -S example -P silico"
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
ALIGN_ARG="-k 14"
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
      if [ ! -z "$OPTARG" ]
      then
        # Use user-input 
        ALIGN_ARG=$OPTARG
      fi
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

# Create Read Group ID
RG="@RG\tLB:$RGLB\tID:$RGID\tSM:$RGSM\tPO:$RGPO\tPL:$RGPL"

# ALIGN ===================================================
# Generate random alpha-numeric for run-id
# RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )

# Logging to STDOUT
echo " -- bowtie2 Alignment Pipeline -- "
echo " date:    $(date)"
echo " version: $PIPE_VERSION "
echo " output:  $OUTNAME"

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
echo " bt2 arguments:   $ALIGN_ARG -p $THREADS"
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
  # -G 12 excludes reads with both reads in pair unmapped
  echo "bwa mem $ALIGN_ARG -t $THREADS \\"
  echo "  -R $RG \\"
  echo "  $GENOME \\"
  echo "  $FQ1 $FQ2 | \\"
  echo "samtools view -bS -G 12 - > "$OUTNAME".bam"

  bwa mem $ALIGN_ARG -t $THREAD \
    -R $RG \
    $GENOME \
    $FQ1 $FQ2 | \
    samtools view -bS -G 12 - > "$OUTNAME".bam

  echo ""
  echo "Alignment complete."

  if [ $DEBUG = "F" ]
  then
    echo "clearing fq-files"
    # Clean-up FQ files to save space
    rm $FQ1 $FQ2
  fi
  # OUTPUT: $OUTNAME.bam
  
else
  # Unpaired Read Alignment
  # -F 4 will retain only mapped reads
  echo "bwa mem $ALIGN_ARG -t $THREADS \\"
  echo "  -R $RG \\"
  echo "  $GENOME \\"
  echo "  $FQ1 $FQ2 | \\"
  echo "samtools view -bS -G 12 - > "$OUTNAME".bam"

  bwa mem $ALIGN_ARG -t $THREAD \
    -R $RG \
    $GENOME \
    $FQ3 | \
    samtools view -bS -F 4 - > "$OUTNAME".bam

  echo ""
  echo "Alignment complete."

    if [ $DEBUG = "F" ]
    then
      echo "clearing fq-files"
      # Clean-up FQ
      rm $FQ3
    fi
  # OUTPUT: $OUTNAME.bam
fi
