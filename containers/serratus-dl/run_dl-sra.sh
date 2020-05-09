#!/usr/bin/bash
set -eu
# Max number of FQ reads to download per accession [100M reads]
# ~ 20 GB max
FQMAX=${FQMAX:-100000000}
# 4*Max number of FQ reads to partition per fq-block [250K reads * 4 lines]
BLOCKSIZE=${BLOCKSIZE:-4000000}
BASEDIR=${BASEDIR:-.}

while getopts n: FLAG; do
  case $FLAG in
    n)
      BLOCKSIZE=$OPTARG
      ;;
    \?) #unrecognized option - show help
      echo "Input parameter not recognized"
      usage
      ;;
  esac
done
shift $((OPTIND-1))

SRA="$1"
S3_BUCKET="$2"

# Outputs created by fastq-dump
FQ_1="$SRA"_1.fastq
FQ_2="$SRA"_2.fastq
FQ_3="$SRA".fastq

# Prefetch the data before processing
# This should be VERY fast, and will cause fastq-dump to have
# smoother CPU usage in the end.
prefetch $SRA

# Create some named pipes for fastq-dump to put its data into.
mkfifo "$FQ_1" "$FQ_2" "$FQ_3"
fastq-dump -X $FQMAX --split-e $SRA & pid=$!

S3_FQ1="s3://$S3_BUCKET/fq-blocks/$SRA/$SRA.1.fq.%010d"
S3_FQ2="s3://$S3_BUCKET/fq-blocks/$SRA/$SRA.2.fq.%010d"
S3_FQ3="s3://$S3_BUCKET/fq-blocks/$SRA/$SRA.3.fq.%010d"

# Stream those pipes into a parallel AWS upload.
# TODO: there are 3 levels of concurrency here: with bash, parallel, and aws.
#       should we disable concurrency in parallel or aws?
cat "$FQ_1" | parallel --block 100M --pipe -N"$BLOCKSIZE" "$BASEDIR"/s3_cp_formatted.sh "$S3_FQ1" "{#}" &
cat "$FQ_2" | parallel --block 100M --pipe -N"$BLOCKSIZE" "$BASEDIR"/s3_cp_formatted.sh "$S3_FQ2" "{#}" &
cat "$FQ_3" | parallel --block 100M --pipe -N"$BLOCKSIZE" "$BASEDIR"/s3_cp_formatted.sh "$S3_FQ3" "{#}" &

# Wait for fastq-dump to finish.
wait $pid

# Send nothing through (open and close the pipes), so that all readers get
# EOF.  Some will still be waiting at this point since fastq-dump doesn't
# open output files unless it has something to write.
echo -n > "$FQ_1"
echo -n > "$FQ_2"
echo -n > "$FQ_3"

# Wait for cat, parallel, etc to get the message.
wait
rm -f "$FQ_1" "$FQ_2" "$FQ_3"

