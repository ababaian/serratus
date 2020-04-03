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
  echo "          Output appended with 'r' suffix"
  echo "    -d    Working directory [$PWD]"
  echo "    -o    <output_filename_prefix>"
  echo "          note: if REVERSE sequence is generated"
  echo "          the suffix 'r' is added to <output_prefix>"
  echo ""
  echo "    Outputs: "
  echo ""
  echo "ex: ./run_mblast-index.sh -x seq/cov0.fa -o cov0 "
  echo "out:      cov0r.fa"
  echo "          cov0r.blastdb.log"
  echo "          cov0r.nhr"
  echo "          cov0r.nin"
  echo "          cov0r.nsq"
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
  TMP="$OUTPUT"r
  OUTNAME=$TMP
else
  TMP="$OUTPUT"
  OUTNAME=$TMP
fi

# ALIGN ===================================================
# Generate random alpha-numeric for run-id
RUNID=$(cat /dev/urandom | tr -dc 'a-z0-9' | fold -w 8 | head -n 1 )
mkdir -p $RUNID; cd $RUNID

# Logging to STDOUT
echo " -- mblast-index Pipeline -- "
echo " date:    $(date)"
echo " version: $PIPE_VERSION "
echo " ami:     $AMI_VERSION  "
echo " output:  $OUTNAME"
echo " run-id:  $RUNID"

# ---------------------------

echo " genome:  $GENOME"
echo " index args : $BT2_ARG -p $THREADS"
echo""
echo 'Initializing ...'
echo ""

if [ "$ADD_REVERSE" = "true" ]
then
  # A control sequence of reversed, non-compliment CoV
  # sequences as internal 'control' for BLAST
  grep -e '^>' $GENOME | tac - | sed 's/>/>REVERSE_/g' - > header.tmp
  sed 's/^>.*/\>/g' $GENOME | tr '\n' '*' | rev | sed 's/>/\n/g' - | \
  paste -d'\n' header.tmp - | sed 's/*/\n/g' - | sed '/^$/d' - > rev.fa.tmp
  rm header.tmp

  # Concatenate reverse sequences into genome file
  GENOMEr=$(echo $GENOME | sed 's/.fa/r.fa/g' -)

  cat $GENOME rev.fa.tmp > $GENOMEr
  rm rev.fa.tmp
  GENOME=$GENOMEr
fi

#TODO layer in taxonomy data to BLAST DB
makeblastdb -dbtype nucl -in $GENOME -title $OUTNAME -out $OUTNAME \
  &> "$OUTNAME".blastdb.log

# output:
#  $OUTNAME.blastdb.log
#  $OUTNAME.nhr
#  $OUTNAME.nin
#  $OUTNAME.nsq
