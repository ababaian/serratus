#!/bin/bash

S3_BASE=https://serratus-public.s3.amazonaws.com

mkdir -p reference
# get the reference alignment and tree, and the 
wget -P reference 	${S3_BASE}/rce/uniprot_genes/pol_msas/pol.muscle.phys \
					${S3_BASE}/rce/uniprot_genes/raxml/pol.muscle/RAxML_bestTree.pol.muscle.raxml \
					${S3_BASE}/rce/complete_cov_genomes/complete.tsv

mkdir -p raw
# get the file specifying which contigs to take
wget -P raw/ ${S3_BASE}/assemblies/analysis/catA-v1.txt 
 
# get the filenames of all cat-A contigs
# and merge the sequences into one fasta file
(while IFS= read -r line; do echo "contigs/${line##*/}"; done < raw/catA-v1.txt) | xargs msa-merge > raw/catA-contigs.fa

# get orfs / individual genes
getorf -sequence raw/catA-contigs.fa -snucleotide1 -sformat1 fasta -outseq raw/orfs.fa -osformat2 fasta
# normalize the orf seq names
sed -i -e "s/[[:space:]]/_/g" raw/orfs.fa

# build the hmm
mkdir -p align
hmmbuild --amino align/ref.hmm reference/pol.muscle.phys

# search orfs against the hmm to get evalues
hmmsearch -o align/search.log --noali --cpu 4 --tblout align/hits.tsv align/ref.hmm raw/orfs.fa

# keep only good hits from the orf file 
seqtk subseq raw/orfs.fa <(grep -v '^#' align/hits.tsv | awk '{print $1}') > align/orfs.filtered.fa

# align the good hits
hmmalign --outformat afa --mapali reference/pol.muscle.phys align/ref.hmm align/orfs.filtered.fa | gzip --best > align/aligned.orfs.afa.gz

# split for epa
mkdir -p place 
epa-ng --outdir place/ --redo --split reference/pol.muscle.phys align/aligned.orfs.afa.gz
gzip --force --best place/query.fasta

# re-get the model params
mkdir -p eval
raxml-ng --evaluate --model PROTGTR+F --tree reference/RAxML_bestTree.pol.muscle.raxml --msa reference/pol.muscle.phys --prefix eval/eval --threads 2

# place
epa-ng --threads 4 --query place/query.fasta.gz --msa place/reference.fasta \
--outdir place/ --model eval/eval.raxml.bestModel --tree eval/eval.raxml.bestTree --redo

mkdir -p assign
# get reference taxonomy file in the right order for gappa assign
# this also fixes the screwed up taxa names to be the same as with the tree (thanks phylip!)
awk -F '\t' '{if(length($1) == 8){$1=sprintf("%s.",$1)};print $1,$6}' OFS='\t' reference/complete.tsv > assign/taxonomy.tsv

# do the assignment
gappa examine assign --jplace-path place/epa_result.jplace --taxon-file assign/taxonomy.tsv --out-dir assign/ --krona --best-hit --per-query-results --allow-file-overwriting

# make per-query best hit results more readable
awk '{split($1,a,".");$1=a[1];print}' OFS='\t'  assign/assign_per_query.tsv > assign/readable.per_query.tsv
