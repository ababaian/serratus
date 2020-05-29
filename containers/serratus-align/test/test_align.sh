#!/usr/bin/env bash
set -euo pipefail

aws s3 cp s3://serratus-public/test-data/fq-blocks/SRR11166696.1.fq.000000.gz ./
aws s3 cp s3://serratus-public/test-data/fq-blocks/SRR11166696.2.fq.000000.gz ./

FQ1='SRR11166696.1.fq.000000.gz'
FQ2='SRR11166696.2.fq.000000.gz'

aws s3 sync s3://serratus-public/seq/cov2m/ ./
GENOME='cov2m'

SRA='UNIT_TEST'
BL_N='1337'

ALIGN_ARGS='--very-sensitive-local'

RGLB='SRAXXX'
RGID='SRBXXX'
RGSM='BIOSXXX'
RGPO='ILLUMINA'

./run_bowtie2.sh \
 -1 $FQ1 -2 $FQ2 -x $GENOME \
 -o $SRA.$BL_N -p 1 -a "$ALIGN_ARGS" \
 -L "$RGLB" -I "$RGID" -S "$RGSM" -P "ILLUMINA"
