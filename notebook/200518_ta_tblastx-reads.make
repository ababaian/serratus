

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


### Utils:

stage-on-s3:
	gzip /media/storage/tblastx/blastx-out.tsv
	aws s3 cp /media/storage/tblastx/blastx-out.tsv.gz s3://serratus-public/notebook/200518_ta_blastx-reads/ERR2756788-blastx-out.tsv.gz
	aws s3 cp /media/storage/tblastx/cov-hit-prots.fsa s3://serratus-public/notebook/200518_ta_blastx-reads/ERR2756788-cov-hit-prots.fsa
