#!/bin/bash

set -e

S3_BASE=https://serratus-public.s3.amazonaws.com
SERRAPLACE=${S3_BASE}/pb/serraplace

OPTIND=1

verbose=0
threads=4

die () {
    echo >&2 "$@"
    exit 1
}

show_help() {
  echo "Usage: $0 [OPTION]... contig_files..."
  echo "Options:"
  printf "  %s\t%s\n" "-h" "show help"
  printf "  %s\t%s\n" "-v" "increase verbosity"
  printf "  %s\t%s\n" "-t" "number of threads"
  printf "  %s\t%s\n" "-c" "alternative catX-file"
}

while getopts "h?vt:c:" opt; do
  case "$opt" in
  h|\?)
    show_help
    exit 0
    ;;
  v)  verbose=1
    ;;
  t)  threads=$OPTARG
    ;;
  c)  catfile=$OPTARG
    ;;
  esac
done

# input validation

# ensure threads is a number
int_regex='^[0-9]+$'
[[ $threads =~ $int_regex ]] || die "Invalid number of threads: $threads"

# ensure catfile exists if it was specified
[[ ! -z $catfile ]] && [[ ! -f "${catfile}" ]] && die "No such file: $catfile"


shift $((OPTIND-1))

wget_mod () {
  [[ $verbose -eq 1 ]] && echo "Ensuring update from $2 to $1"
  wget -qNO $1 $2
}

REFPHY=reference/pol.reduced.phy
TREE=reference/pol.reduced.newick
MODEL=reference/pol.reduced.raxml.bestModel
TAXONOMY=reference/complete.tsv

mkdir -p reference
# get the reference alignment, model and tree, and the taxonomy file
wget_mod ${REFPHY} ${SERRAPLACE}/tree/pol.muscle.phys.reduced
wget_mod ${MODEL} ${SERRAPLACE}/tree/pol.reduced.raxml.bestModel
wget_mod ${TREE} ${SERRAPLACE}/tree/pol.reduced.raxml.bestTree
wget_mod ${TAXONOMY} ${S3_BASE}/rce/complete_cov_genomes/complete.tsv

mkdir -p raw
# if contig files were not passed via command line, download them from the specified file
if [[ $# -eq 0 ]]
then
  if [[ -z catfile ]]
  then 
    CATX=raw/catX-spec.txt
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
  msa-merge $@ > raw/contigs.fa
fi

# get orfs / individual genes
getorf -sequence raw/contigs.fa -snucleotide1 -sformat1 fasta -outseq raw/orfs.fa -osformat2 fasta
# normalize the orf seq names
sed -i -e "s/[[:space:]]/_/g" raw/orfs.fa

mkdir -p align
# how to build the hmm:
# hmmbuild --amino align/ref.hmm ${REFPHY}
# but we will download it instead
REF_HMM=align/ref.hmm
wget_mod ${REF_HMM} ${SERRAPLACE}/reference/ref.hmm

# search orfs against the hmm to get evalues
hmmsearch -o align/search.log --noali --cpu ${threads} --tblout align/hits.tsv ${REF_HMM} raw/orfs.fa

# keep only good hits from the orf file 
seqtk subseq raw/orfs.fa <(grep -v '^#' align/hits.tsv | awk '{print $1}') > align/orfs.filtered.fa

# align the good hits
hmmalign --outformat afa --mapali ${REFPHY} ${REF_HMM} align/orfs.filtered.fa | gzip --best > align/aligned.orfs.afa.gz

# split for epa
mkdir -p place 
epa-ng --outdir place/ --redo --split ${REFPHY} align/aligned.orfs.afa.gz
gzip --force --best place/query.fasta

# place
epa-ng --threads ${threads} --query place/query.fasta.gz --msa place/reference.fasta \
--outdir place/ --model ${MODEL} --tree ${TREE} --redo --no-heur

mkdir -p assign
# get reference taxonomy file in the right order for gappa assign
# this also fixes the screwed up taxa names to be the same as with the tree (thanks phylip!)
awk -F '\t' '{if(length($1) == 8){$1=sprintf("%s.",$1)};print $1,$6}' OFS='\t' ${TAXONOMY} > assign/taxonomy.tsv

# do the assignment
gappa examine assign --jplace-path place/epa_result.jplace --taxon-file assign/taxonomy.tsv \
--out-dir assign/ --per-query-results --allow-file-overwriting --consensus-thresh 0.66 --log-file assign/assign.log

# make per-query best hit results more readable
awk '{split($1,a,".");$1=a[1];print}' OFS='\t'  assign/assign_per_query.tsv > assign/readable.per_query.tsv

# if this is a bigger job, produce a grafted tree
if [[ $# -eq 0 ]]
then
  gappa examine graft --jplace-path place/epa_result.jplace --name-prefix "SERRATUS_" --out-dir assign/ \
  --threads $threads --log-file assign/graft.log --redo
fi
