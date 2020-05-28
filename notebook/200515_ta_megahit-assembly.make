

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
DATA_DIR ?= /media/ephemeral/assembly
OUT_DIR  ?= $(CURDIR)/
INSTALL_DIR ?= $(CURDIR)/third-party

LD_LIBRARY_PATH = /usr/local/lib
export LD_LIBRARY_PATH

PATH := $(PATH):$(INSTALL_DIR)/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$(INSTALL_DIR)/ncbi-blast-2.10.0+/bin:$(INSTALL_DIR)/edirect
export PATH

## Previously-assembled samples:
# ERR2756787 \
# 	ERR3569452 \
# 	ERR2756788 \
# 	SRR7287114 \
# 	SRR7287110 \
# 	SRR10829950 \

SRA-list := \
	SRR924370 \
	SRR11616465

SRAs := $(shell echo $(SRA-list) | tr ' ' '\n')

SRA-flags  := $(addsuffix .txt, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-fastqs := $(addsuffix _1.fastq, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-assemblies := $(addsuffix _megahit/final.contigs.fa, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-alignments := $(addsuffix .delta, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-align-reports:= $(addsuffix .delta, $(addprefix $(DATA_DIR)/, $(SRAs)))

#### Targets:

### Debug:
print-env:
	echo $(SRAs)
	echo $(SRA-fastqs)
	echo $(SRA-flags)
	echo $(notdir $(basename $(SRA-flags)))
	echo $(SRA-assemblies)
	echo $(SRA-alignments)

### Installation:
install-ubuntu-deps:
	sudo apt-get update
	sudo apt-get install --assume-yes awscli libcurl4-openssl-dev cpanminus perl-doc parallel libssl-dev zlib1g-dev

install-megahit:
	mkdir -p $(INSTALL_DIR)
	cd $(INSTALL_DIR)
	wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
	tar xzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

test-megahit:
	cd $(INSTALL_DIR)/MEGAHIT-1.2.9-Linux-x86_64-static
	bin/megahit --test

$(INSTALL_DIR)/mummer:
	mkdir -p $(INSTALL_DIR)
	cd $(INSTALL_DIR)
	wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
	tar xzf mummer-4.0.0beta2.tar.gz
	cd mummer-4.0.0beta2
	./configure
	make
	sudo make install

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


$(INSTALL_DIR)/prokka:
	mkdir -p $(INSTALL_DIR)
	sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
	sudo cpan Bio::Perl
	cd $(INSTALL_DIR)
	git clone https://github.com/tseemann/prokka.git
	prokka/bin/prokka --setupdb

$(INSTALL_DIR)/eggnog-mapper:
	mkdir -p $(INSTALL_DIR)
	cd $(INSTALL_DIR)
	git clone git@github.com:eggnogdb/eggnog-mapper.git

### Stage Data:
setup-runs:
	mkdir -p $(DATA_DIR)
	touch $(SRA-flags)

stage-data: $(DATA_DIR)/cov2m.fa $(HOME)/.ncbirc $(SRA-fastqs)

$(DATA_DIR)/%_1.fastq: $(DATA_DIR)/%.txt
	mkdir -p $(DATA_DIR)	
	cd $(DATA_DIR)
	fastq-dump -I --split-files --skip-technical $(basename $(notdir $^))
	if [ -e $(basename $(notdir $^))_3.fastq ]
	then
		mv $(basename $(notdir $^))_3.fastq $(basename $(notdir $^))_1.fastq
	fi
	if [ -e $(basename $(notdir $^))_4.fastq ]
	then
		mv $(basename $(notdir $^))_4.fastq $(basename $(notdir $^))_2.fastq
	fi

# $(DATA_DIR)/ERR2756788_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files ERR2756788

# $(DATA_DIR)/ERR2756787_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files ERR2756787

# $(DATA_DIR)/ERR3569452_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files ERR3569452

# $(DATA_DIR)/SRR7287114_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files SRR7287114

# $(DATA_DIR)/SRR7287110_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files SRR7287110

# $(DATA_DIR)/SRR10829950_1.fastq:
# 	mkdir -p $(DATA_DIR)
# 	cd $(DATA_DIR)
# 	fastq-dump -I --split-files SRR10829950

$(DATA_DIR)/cov2m.fa:
	aws s3 cp s3://serratus-public/seq/cov2m/cov2m.fa $(DATA_DIR)/cov2m.fa

$(HOME)/.ncbirc:
	@echo "[BLAST]" > $(HOME)/.ncbirc
	@echo "BLASTDB=$(DATA_DIR)/../blastdbs" >> $(HOME)/.ncbirc

update-blastdbs:
	mkdir -p $(DATA_DIR)/../blastdbs
	cd $(DATA_DIR)/../blastdbs
	date; time update_blastdb.pl --decompress nt; echo $$?; date

download-emapper-db:
	$(INSTALL_DIR)/eggnog-mapper/download_eggnog_data.py \
		-y -f --data_dir $(DATA_DIR)

### Run Assembly:
## ~30 minutes on 72 cores, 144G RAM, and 200G instance swap
run-assembly: $(SRA-assemblies)

$(DATA_DIR)/%_megahit/final.contigs.fa: $(DATA_DIR)/%_1.fastq
	if [ -e `echo $^ | sed 's/_1.fastq/_2.fastq/'` ]
	then	
		megahit -1 $^ \
			-2 `echo $^ | sed 's/_1.fastq/_2.fastq/' ` \
			-o `echo $^ | sed 's/_1.fastq/_megahit/' ` \
			--mem-flag 2
	else
		megahit -r $^ \
			-o `echo $^ | sed 's/_1.fastq/_megahit/' ` \
			--mem-flag 2
	fi

# $(DATA_DIR)/ERR2756788_megahit/final.contigs.fa: $(DATA_DIR)/ERR2756788_1.fastq
# 	megahit -1 $(DATA_DIR)/ERR2756788_1.fastq \
# 		-2 $(DATA_DIR)/ERR2756788_2.fastq \
# 		-o $(DATA_DIR)/megahit_out \
# 		--mem-flag 2

# $(DATA_DIR)/ERR2756787_megahit/final.contigs.fa: $(DATA_DIR)/ERR2756787_1.fastq
# 	megahit -1 $(DATA_DIR)/ERR2756787_1.fastq \
# 		-2 $(DATA_DIR)/ERR2756787_2.fastq \
# 		-o $(DATA_DIR)/ERR2756787_megahit \
# 		--mem-flag 2

# $(DATA_DIR)/ERR3569452_megahit/final.contigs.fa: $(DATA_DIR)/ERR3569452_1.fastq
# 	megahit -1 $(DATA_DIR)/ERR3569452_1.fastq \
# 		-2 $(DATA_DIR)/ERR3569452_2.fastq \
# 		-o $(DATA_DIR)/ERR3569452_megahit \
# 		--mem-flag 2

# $(DATA_DIR)/SRR7287114_megahit/final.contigs.fa: $(DATA_DIR)/SRR7287114_1.fastq
# 	megahit -1 $(DATA_DIR)/SRR7287114_1.fastq \
# 		-2 $(DATA_DIR)/SRR7287114_2.fastq \
# 		-o $(DATA_DIR)/SRR7287114_megahit \
# 		--mem-flag 2

# $(DATA_DIR)/SRR7287110_megahit/final.contigs.fa: $(DATA_DIR)/SRR7287110_1.fastq
# 	megahit -1 $(DATA_DIR)/SRR7287110_1.fastq \
# 		-2 $(DATA_DIR)/SRR7287110_2.fastq \
# 		-o $(DATA_DIR)/SRR7287110_megahit \
# 		--mem-flag 2

# $(DATA_DIR)/SRR10829950_megahit/final.contigs.fa: $(DATA_DIR)/SRR10829950_1.fastq
# 	megahit -1 $(DATA_DIR)/SRR10829950_1.fastq \
# 		-2 $(DATA_DIR)/SRR10829950_2.fastq \
# 		-o $(DATA_DIR)/SRR10829950_megahit \
# 		--mem-flag 2

### Assembly Alignment:

align-assembly: $(SRA-alignments)

$(DATA_DIR)/%.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/%_megahit/final.contigs.fa
	cd $(DATA_DIR)
	nucmer --delta=$@ $^

# $(DATA_DIR)/ERR3569452.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR3569452_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=ERR3569452 $^

# $(DATA_DIR)/ERR2756787.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR2756787_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=ERR2756787 $^

# $(DATA_DIR)/ERR2756788.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR2756788_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=ERR2756788 $^

# $(DATA_DIR)/SRR7287114.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/SRR7287114_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=SRR7287114 $^

# $(DATA_DIR)/SRR7287110.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/SRR7287110_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=SRR7287110 $^

# $(DATA_DIR)/SRR10829950.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/SRR10829950_megahit/final.contigs.fa
# 	cd $(DATA_DIR)
# 	nucmer --prefix=SRR10829950 $^


## This shows that k141_146741 wholly contains three fragments in cov2m, and is ~29k:
show-coords:
	cd $(DATA_DIR)
	for alignment in *.delta
	do
		show-coords -c -q -g -o -l -T \
			$$alignment \
			> `echo $$alignment | sed 's/.delta/-alignments.txt/'`
	done

## Blast:

$(DATA_DIR)/%/blast.tsv: $(DATA_DIR)/%_megahit/final.contigs.fa
	cd $(DATA_DIR)
	mkdir -p `dirname $@`
	get_species_taxids.sh -t 10239 > virus-ids.txt
	blastn -query $^ \
		-db nt \
		-evalue 1 \
		-task blastn \
		-dust yes \
		-taxidlist virus-ids.txt \
		-num_threads `nproc` \
		-outfmt 7 \
		> $@

run-blast: $(DATA_DIR)/ERR2756788/blast.tsv




### Utils:

stage-results-to-s3:
	cd $(DATA_DIR)
	for acc in $(SRAs)
	do
		[ ! -e $${acc}_megahit.tar.gz ] && tar czf $${acc}_megahit.tar.gz $${acc}_megahit
		aws s3 cp $${acc}_megahit.tar.gz s3://serratus-public/notebook/200515_ta/
	done

stage-alignment-reports:
	cp $(DATA_DIR)/*-alignments.txt $(CURDIR)/notebook/200515_ta_megahit-assembly/

## Use lsblk to find a suitable ephemeral drive to use, then call it like:
## SWAP_DEVICE=/dev/nvme1n1 SWAP_SIZE=400G make set-up-swap
set-up-swap:
	sudo mkdir -p /media/ephemeral
	sudo mkfs.ext4 $(SWAP_DEVICE)
	sudo mount $(SWAP_DEVICE) /media/ephemeral
	sudo chmod -R 777 /media/ephemeral
	sudo fallocate -l $(SWAP_SIZE) /media/ephemeral/swapfile
	sudo chmod 600 /media/ephemeral/swapfile
	sudo mkswap /media/ephemeral/swapfile
	sudo swapon /media/ephemeral/swapfile
