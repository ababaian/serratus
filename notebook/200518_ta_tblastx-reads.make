

#### Assumptions:
## * sudo apt-get install make
## * AWS CLI configured


#### Make configuration:

## Use Bash as default shell, and in strict mode:
SHELL := /bin/bash
.SHELLFLAGS = -ec

## If the parent env doesn't ste TMPDIR, do it ourselves:
TMPDIR ?= /tmp

## This makes all recipe lines execute within a shared shell process:
## https://www.gnu.org/software/make/manual/html_node/One-Shell.html#One-Shell
.ONESHELL:

## If a recipe contains an error, delete the target:
## https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#Special-Targets
.DELETE_ON_ERROR:

## This is necessary to make sure that these intermediate files aren't clobbered:
.SECONDARY:

#### Local Parameters:
DATA_DIR ?= /media/storage/tblastx
OUT_DIR  ?= $(CURDIR)/
INSTALL_DIR ?= $(CURDIR)/third-party

PATH := $(PATH):$(INSTALL_DIR)/edirect:$(INSTALL_DIR)/ncbi-blast-2.10.0+/bin
export PATH

#### Targets:

### Installation:

install-ubuntu-deps:
	sudo apt-get update
	sudo apt-get install --assume-yes awscli libcurl4-openssl-dev cpanminus perl-doc parallel libssl-dev zlib1g-dev

$(INSTALL_DIR)/edirect:
	cd $(INSTALL_DIR)
	sudo cpanm IO::Socket::SSL
	sudo cpanm LWP::Protocol::https
	perl -MNet::FTP -e \
		'$$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
		$$ftp->login; $$ftp->binary;
		$$ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
	gunzip -c edirect.tar.gz | tar xf -
	rm edirect.tar.gz
	./edirect/setup.sh

$(INSTALL_DIR)/ncbi-blast-2.10.0+/bin/blastn: third-party/edirect
	mkdir -p $(INSTALL_DIR)
	cd $(INSTALL_DIR)
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz
	tar xzf ncbi-blast-2.10.0+-x64-linux.tar.gz


### Data Prep:
$(DATA_DIR)/cov-hit-prots.fsa:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	aws s3 cp s3://serratus-public/out/200505_zoonotic/summary/ERR2756788.summary .
	awk -F';' '/^acc=/ {split($$1,parts,"="); if(parts[2] != "pan_genome") print parts[2]}' \
		ERR2756788.summary \
	| epost -db nuccore \
		-label coronaprot \
		| elink -target protein \
		| efetch -format fasta \
		> cov-hit-prots.fsa

$(DATA_DIR)/cov-hit-prots.fsa.pdb:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	makeblastdb -in cov-hit-prots.fsa \
		-parse_seqids \
		-blastdb_version 5 \
		-title "CoV Hit Proteins" \
		-dbtype prot

$(DATA_DIR)/query-reads.fa:
	awk '(NR%4)==1 { gsub(/^@/,">"); print; next } \
             (NR%4)==2' \
		$(DATA_DIR)/../assembly/ERR2756788_1.fastq \
		$(DATA_DIR)/../assembly/ERR2756788_2.fastq \
		> $@

### Run:

parallel-blastx: $(DATA_DIR)/query-reads.fa
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	cat $^ \
		| parallel --pipe --progress -N1000 --recstart '>' blastx -db $(DATA_DIR)/cov-hit-prots.fsa -evalue .001 -outfmt 7 \
		> blastx-out.tsv

run-blastx: $(DATA_DIR)/query-reads.fa
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	blastx -query $^ \
		-db $(DATA_DIR)/cov-hit-prots.fsa \
		-evalue .001 \
		-num_threads `nproc` \
		-outfmt 7 \
		> $@

## megablast: 2m26s
## blastn: 

eval-blastn_vdb-mapping:
	cd /dev/shm/
	mkdir -p blastn-map
	cd blastn-map	
	wget https://github.com/ababaian/serratus/wiki/assets/Fr4NK.fa
	revseq -sequence Fr4NK.fa -outseq Fr4NK_revcomp.fa
	samtools faidx Fr4NK_revcomp.fa
	#wget -r https://sra-download.ncbi.nlm.nih.gov/traces/era23/ERR/ERR2756/ERR2756788
	#cd sra-download.ncbi.nlm.nih.gov/traces/era23/ERR/ERR2756
	time blastn_vdb \
		-db "ERR2756788" \
		-query Fr4NK_revcomp.fa \
		-evalue .001 \
		-num_threads `nproc` \
		-outfmt 17 \
		-task megablast \
		-max_target_seqs 100000 \
		| sed 's/Query_1/Bat/' \
		| samtools view -S -b \
		| samtools sort \
		> test-megablast-sorted.bam
	samtools index test-megablast-sorted.bam
	time blastn_vdb \
		-db "ERR2756788" \
		-query Fr4NK_revcomp.fa \
		-evalue .001 \
		-num_threads `nproc` \
		-outfmt 17 \
		-task blastn \
		-max_target_seqs 100000 \
		| sed 's/Query_1/Bat/' \
		| samtools view -S -b \
		| samtools sort \
		> test-blastn-sorted.bam

## using default -> 7892 hard clipping
## using 5/-4/8/6 -> 7323 hard clipping
## using 1/-1/0/2 -> 
eval-blastn_vdb-variant-calling:
	cd /dev/shm/blastn-map
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/503/155/GCF_001503155.1_ViralProj307859/GCF_001503155.1_ViralProj307859_genomic.fna.gz
	gunzip GCF_001503155.1_ViralProj307859_genomic.fna.gz
	prefetch -p "ERR2756788"
	time blastn_vdb \
		-db "ERR2756788" \
		-query Fr4NK_revcomp.fa \
		-evalue 10 \
		-num_threads `nproc` \
		-outfmt "17 SQ" \
		-task blastn \
		-gapopen 8 \
		-gapextend 6 \
		-penalty -4 \
		-reward 5 \
		-max_target_seqs 100000 \
		| sed 's/Query_1/NC_028811.1/' \
		| samtools view -S -b \
		| samtools sort \
		> frank-variants-sorted.bam
	samtools index frank-variants-sorted.bam

## Need to find hard clipped alignments, record the hits, the hard-clip numbers,
## and then search the FASTQ file to find the sequence & 
soft-clip-to-hard-clip:
	samtools view frank-variants-sorted.bam \
		| awk -F"\t" '{split($$1, vdb_coord, "."); print vdb_coord[2]}' \
		| sort | uniq \
		> spots.txt
	vdb-dump -f tab -C SPOT_ID,NAME,READ_START,READ,QUALITY
	export LC_ALL=C
	awk -F"\t" -v OFS="\t" \
		'{ split($$1, vdb_coord, ".");
		   cmd = "vdb-dump -f tab -C SPOT_ID,NAME,READ_START,READ,QUALITY -R " vdb_coord[2] " ERR2756788";
		   cmd | getline vdb_line;
		   split(vdb_line, vdb_stats);
		   split(vdb_stats[3],start_coords,", ");
		   num_scores = split(vdb_stats[5],q_scores,", ");
		   q_str = "";		
		   if( vdb_coord[3] == "1" ) {
			$10 = substr(vdb_stats[4],1,start_coords[2]);
			for(i=start_coords[1]+1;i<=start_coords[2];i++)
				q_str = q_str sprintf("%c", 33+q_scores[i]);
			$$11 = q_str;
		   }
		   else {
			$$10 = substr(vdb_stats[4],start_coords[2]+1);
			for(i=start_coords[2]+1;i<=length(q_scores);i++)
				q_str = q_str sprintf("%c", 33+q_scores[i]);
			$$11 = q_str;
		   }
		   print;
		 }' \
	<(samtools view frank-variants-sorted.bam)

bowtie2-align:
	cd /dev/shm/blastn-map
	prefetch -p "ERR2756788"
	cd /tmp
	fastq-dump --split-e "ERR2756788"
	cd /dev/shm/blastn-map
	time bowtie2 --very-sensitive-local -x /dev/shm/blastn-map/NC_028811 -U /tmp/ERR2756788_1.fastq,/tmp/ERR2756788_2.fastq -S > /dev/shm/blastn-map/bt2-align-very-sensitive-local.sam

### Utils:

stage-on-s3:
	gzip /media/storage/tblastx/blastx-out.tsv
	aws s3 cp /media/storage/tblastx/blastx-out.tsv.gz s3://serratus-public/notebook/200518_ta_blastx-reads/ERR2756788-blastx-out.tsv.gz
	aws s3 cp /media/storage/tblastx/cov-hit-prots.fsa s3://serratus-public/notebook/200518_ta_blastx-reads/ERR2756788-cov-hit-prots.fsa












































