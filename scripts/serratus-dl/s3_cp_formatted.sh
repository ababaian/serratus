#!/usr/bin/bash
set -eu

# Run `aws s3 cp`, but format the output using printf first, and read the
# input file from a pipe
#
#     $ cat my_file | s3_ul_formatted s3://bucket/sra_%04d.fastq 42
#   -> uploads my_file to the an object called sra_0042.fastq to S3.
#
# Args 1 and 2 will be passed through printf.  All other args
# will be passed to `aws s3 cp`
#
# This script is intended to be an argument to the parallel command, like:
#
#     $ cat my_pipe | parallel --pipe -N100 s3_cp_formatted s3://bucket/file_%010d.fastq {#}
# 
# Which will create:
#     file_0000000000.fastq
#     file_0000000001.fastq
#     [...]
# on S3, in parallel.
FMT="$1"
N=$(expr "$2" - 1) # Parallel is one based. :/

DEST=$(printf "$FMT" "$N")
shift 2
aws s3 cp - "$DEST" "$@" <&0
