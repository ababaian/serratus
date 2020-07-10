#!/bin/bash

set -e

S3_BASE=https://serratus-public.s3.amazonaws.com
SERRAPLACE=${S3_BASE}/pb/serraplace

OPTIND=1

verbose=0
threads=4
no_merge=0
download=0

die () {
    echo >&2 "ABORT: $@"
    exit 1
}

show_help() {
  echo "Usage: $0 [OPTION]... contig_files..."
  echo "Options:"
  printf "  %s\t%s\n" "-h" "show help"
  printf "  %s\t%s\n" "-v" "increase verbosity"
  printf "  %s\t%s\n" "-t" "number of threads"
  printf "  %s\t%s\n" "-d" "re-download reference data"
  printf "  %s\t%s\n" "-c" "alternative catX-file"
  printf "  %s\t%s\n" "-m" "turn off merging of explicitly passed contig file (assumes one file with sensical fasta names)"
}

while getopts "h?vmt:c:d" opt; do
  case "$opt" in
  h|\?)
    show_help
    exit 0
    ;;
  v)  verbose=1
    ;;
  d)  download=1
    ;;
  m)  no_merge=1
    ;;
  t)  threads=$OPTARG
    ;;
  c)  catfile=$OPTARG
    ;;
  esac
done
shift $((OPTIND-1))

# input validation

# ensure threads is a number
int_regex='^[0-9]+$'
[[ $threads =~ $int_regex ]] || die "Invalid number of threads: $threads"

# ensure catfile exists if it was specified
[[ ! -z $catfile ]] && [[ ! -f "${catfile}" ]] && die "No such file: $catfile"

[[ $no_merge -eq 1 ]] && [[ $# -ne 1 ]] && die "Turned off merging but specified more than one (or no) file"

wget_mod () {
  [[ $verbose -eq 1 ]] && echo "Ensuring update from $2 to $1"
  wget -qNO $1 $2
}

expect_file () {
  if [[ ! -f $1 ]]
  then
    die "File missing: $1"
  fi
}

REF_MSA=reference/reference.afa
TREE=reference/raxml.bestTree
MODEL=reference/raxml.bestModel
TAXONOMY=reference/complete.tsv
OUTGROUP=reference/outgroupspec.txt

if [[ $download -eq 1 ]]
then
  # get the reference alignment, model and tree, and the taxonomy file
  wget_mod ${REF_MSA} ${SERRAPLACE}/reference/tree/clust.comb.afa
  wget_mod ${MODEL} ${SERRAPLACE}/reference/tree/10_search.raxml.bestModel
  wget_mod ${TREE} ${SERRAPLACE}/reference/tree/10_search.raxml.bestTree
  wget_mod ${OUTGROUP} ${SERRAPLACE}/reference/tree/outgroupspec.txt
  wget_mod ${TAXONOMY} ${SERRAPLACE}/reference/complete.tsv
fi

expect_file ${REF_MSA}
expect_file ${MODEL}
expect_file ${TREE}
expect_file ${OUTGROUP}
expect_file ${TAXONOMY}


mkdir -p reference
# get the reference alignment, model and tree, and the taxonomy file
# wget_mod ${REF_MSA} ${SERRAPLACE}/reference/tree/clust.comb.afa
# wget_mod ${MODEL} ${SERRAPLACE}/reference/tree/10_search.raxml.bestModel
# wget_mod ${TREE} ${SERRAPLACE}/reference/tree/10_search.raxml.bestTree
# wget_mod ${OUTGROUP} ${SERRAPLACE}/reference/tree/outgroupspec.txt
# wget_mod ${TAXONOMY} ${SERRAPLACE}/reference/complete.tsv


mkdir -p raw
CONTIGS=raw/contigs.fa
# if contig files were not passed via command line, download them from the specified file
if [[ $# -eq 0 ]]
then
  CATX=raw/catX-spec.txt
  if [[ -z $catfile ]]
  then
    # get the file specifying which contigs to take
    wget_mod ${CATX} ${S3_BASE}/assemblies/analysis/catA-v1.txt
  else
    CATX=$catfile
  fi

  # if there already is a contigs folder, use that. else download the files specified in the catX-file
  if [[ ! -d contigs/ ]]
  then
    echo "Downloading contigs since I didn't find a contigs/ folder"
    mkdir contigs
    while IFS= read -r line;
    do
      wget_mod "contigs/${line##*/}" "${S3_BASE}/assemblies/contigs/${line##*/}";
    done < ${CATX}
  fi

  # get the filenames of all cat-A contigs
  # and merge the sequences into one fasta file
  (while IFS= read -r line; do echo "contigs/${line##*/}"; done < ${CATX}) | xargs msa-merge > raw/contigs.fa
# if they were passed, just parse those in
else
  if [[ $no_merge -eq 0 ]] 
  then
    msa-merge $@ > ${CONTIGS}
  else
    CONTIGS=$@
    echo "Selected single combined contig file: ${CONTIGS}"
  fi
fi

# get orfs / individual genes
esl-translate ${CONTIGS} > raw/orfs.fa

# normalize the orf seq names
sed -i -e "s/[[:space:]]/_/g" raw/orfs.fa

mkdir -p align
REF_HMM=align/ref.hmm
# how to build the hmm:
# hmmbuild --amino ${REF_HMM} ${REF_MSA}
# but we will download it instead, if we want to
[[ $download -eq 1 ]] && wget_mod ${REF_HMM} ${SERRAPLACE}/reference/hmm/clust.comb.hmm

expect_file ${REF_HMM}

# search orfs against the hmm to get evalues
echo "Running hmmsearch"
hmmsearch -o align/search.log --noali -E 0.01 --cpu ${threads} --tblout align/hits.tsv ${REF_HMM} raw/orfs.fa

# keep only good hits from the orf file 
seqtk subseq raw/orfs.fa <(grep -v '^#' align/hits.tsv | awk '{print $1}') > raw/orfs.filtered.fa
sed -i -E 's/^>.*source=/>/g' raw/orfs.filtered.fa
sed -i -E 's/\..*=$//g' raw/orfs.filtered.fa

# align the good hits
echo "Running hmmalign"
hmmalign --outformat afa --mapali ${REF_MSA} ${REF_HMM} raw/orfs.filtered.fa | gzip --best > align/aligned.orfs.afa.gz

# split for epa
mkdir -p place 
epa-ng --outdir place/ --redo --split ${REF_MSA} align/aligned.orfs.afa.gz
gzip --force --best place/query.fasta

# place
epa-ng --threads ${threads} --query place/query.fasta.gz --msa place/reference.fasta \
--outdir place/ --model ${MODEL} --tree ${TREE} --redo --no-heur

mkdir -p assign
# get reference taxonomy file in the right order for gappa assign
# this also fixes the screwed up taxa names to be the same as with the tree (thanks phylip!)
awk -F '\t' '{print $1,$6}' OFS='\t' ${TAXONOMY} > assign/taxonomy.tsv

# do the assignment
gappa examine assign --jplace-path place/epa_result.jplace --taxon-file assign/taxonomy.tsv \
--out-dir assign/ --per-query-results --allow-file-overwriting --consensus-thresh 0.66 --log-file assign/assign.log \
--root-outgroup ${OUTGROUP} --threads ${threads}

# make per-query best hit results more readable
awk '{split($1,a,".");$1=a[1];print}' OFS='\t'  assign/assign_per_query.tsv > assign/readable.per_query.tsv

# if this is a bigger job, produce a grafted tree
if [[ $# -eq 0 ]]
then
  gappa examine graft --jplace-path place/epa_result.jplace --name-prefix "SERRATUS_" --out-dir assign/ \
  --threads ${threads} --log-file assign/graft.log --redo
fi
