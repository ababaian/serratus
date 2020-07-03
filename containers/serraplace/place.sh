#!/bin/bash

set -e

S3_BASE=https://serratus-public.s3.amazonaws.com
SERRAPLACE=${S3_BASE}/pb/serraplace

OPTIND=1

while getopts "h?" opt; do
  case "$opt" in
  h|\?)
    echo "Usage: $0 [-h] contig_files..."
    exit 0
    ;;
  # v)  verbose=1
  #   ;;
  esac
done

shift $((OPTIND-1))

mkdir -p reference
# get the reference alignment, model and tree, and the taxonomy file
wget -qP reference  ${SERRAPLACE}/tree/pol.muscle.phys.reduced \
                    ${SERRAPLACE}/tree/pol.reduced.raxml.bestModel \
                    ${SERRAPLACE}/tree/pol.reduced.raxml.bestTree \
                    ${S3_BASE}/rce/complete_cov_genomes/complete.tsv

TREE=reference/pol.reduced.raxml.bestTree
MODEL=reference/pol.reduced.raxml.bestModel
REFPHY=reference/pol.muscle.phys.reduced

mkdir -p raw
# if contig files were not passed via command line, download them from the specified file
if [[ $# -eq 0 ]]
then
  # get the file specifying which contigs to take
  wget -qP raw/ ${S3_BASE}/assemblies/analysis/catA-v1.txt 

  # if there already is a contigs folder, use that. else download the files specified in the catX-file
  if [[ ! -d contigs/ ]]
  then
    echo "Downloading contigs since I didn't find a contigs/ folder"
    mkdir contigs
    while IFS= read -r line;
    do
      wget -qP contigs "${S3_BASE}/assemblies/contigs/${line##*/}";
    done < raw/catA-v1.txt
  fi

  # get the filenames of all cat-A contigs
  # and merge the sequences into one fasta file
  (while IFS= read -r line; do echo "contigs/${line##*/}"; done < raw/catA-v1.txt) | xargs msa-merge > raw/contigs.fa
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
wget -qP align ${SERRAPLACE}/reference/ref.hmm

# search orfs against the hmm to get evalues
hmmsearch -o align/search.log --noali --cpu 4 --tblout align/hits.tsv align/ref.hmm raw/orfs.fa

# keep only good hits from the orf file 
seqtk subseq raw/orfs.fa <(grep -v '^#' align/hits.tsv | awk '{print $1}') > align/orfs.filtered.fa

# align the good hits
hmmalign --outformat afa --mapali ${REFPHY} align/ref.hmm align/orfs.filtered.fa | gzip --best > align/aligned.orfs.afa.gz

# split for epa
mkdir -p place 
epa-ng --outdir place/ --redo --split ${REFPHY} align/aligned.orfs.afa.gz
gzip --force --best place/query.fasta

# place
epa-ng --threads 4 --query place/query.fasta.gz --msa place/reference.fasta \
--outdir place/ --model ${MODEL} --tree ${TREE} --redo --no-heur

mkdir -p assign
# get reference taxonomy file in the right order for gappa assign
# this also fixes the screwed up taxa names to be the same as with the tree (thanks phylip!)
awk -F '\t' '{if(length($1) == 8){$1=sprintf("%s.",$1)};print $1,$6}' OFS='\t' reference/complete.tsv > assign/taxonomy.tsv

# do the assignment
gappa examine assign --jplace-path place/epa_result.jplace --taxon-file assign/taxonomy.tsv \
--out-dir assign/ --per-query-results --allow-file-overwriting --consensus-thresh 0.66

# make per-query best hit results more readable
awk '{split($1,a,".");$1=a[1];print}' OFS='\t'  assign/assign_per_query.tsv > assign/readable.per_query.tsv
