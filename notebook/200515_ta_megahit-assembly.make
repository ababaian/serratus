

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
#	SRR924370 \
#	SRR11616465

SRA_LIST := SRR10951638 SRR10951639 SRR10951640 SRR10951641 SRR10951642 SRR10951643 SRR10951644 SRR10951645 SRR10951646 SRR10951647 SRR10951648 SRR10951649 SRR10951650 SRR10951651 SRR10951652 SRR10951653 SRR10951666 SRR10951667 SRR10951668 SRR10951669 SRR10951670 SRR10951671 SRR10951672 SRR10951673 SRR10951674 SRR10951675 SRR10951676 SRR10951677 SRR10951678 SRR10951679 SRR10951680 SRR10951681 SRR10951682 SRR10951684 SRR10951685 SRR10951686 SRR10951687 SRR10951654 SRR10951655 SRR10951657 SRR10951658 SRR10951659 SRR10951660 SRR10951661 SRR10951662 SRR10951663 SRR10951664 SRR10951665 SRR10951656 SRR10951683

SRAs := $(shell echo $(SRA_LIST) | tr ' ' '\n')

SRA-flags  := $(addsuffix .txt, $(addprefix $(DATA_DIR)/fastq/, $(SRAs)))
SRA-fastqs := $(addsuffix _1.fastq, $(addprefix $(DATA_DIR)/fastq/, $(SRAs)))
SRA-assemblies := $(addsuffix _megahit/final.contigs.fa, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-alignments := $(addsuffix .delta, $(addprefix $(DATA_DIR)/, $(SRAs)))
SRA-align-reports:= $(addsuffix .delta, $(addprefix $(DATA_DIR)/, $(SRAs)))

#### Targets:

### Debug:
print-env:
	echo $(SRA_LIST)
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
	mkdir -p $(DATA_DIR)/fastq
	touch $(SRA-flags)


stage-pig-data:
	mkdir -p $(DATA_DIR)/pigpen
	cd $(DATA_DIR)/pigpen
	prefetch `esearch -db sra -query PRJNA602689 | efetch -format uid`

stage-data: $(DATA_DIR)/cov2m.fa $(HOME)/.ncbirc $(SRA-fastqs)

$(DATA_DIR)/%_1.fastq: $(DATA_DIR)/%.txt
	mkdir -p $(DATA_DIR)/fastq	
	cd $(DATA_DIR)/fastq
	fastq-dump -I --split-files --skip-technical $(basename $(notdir $^))
	if [ -e $(basename $(notdir $^))_3.fastq ]
	then
		mv $(basename $(notdir $^))_3.fastq $(basename $(notdir $^))_1.fastq
	fi
	if [ -e $(basename $(notdir $^))_4.fastq ]
	then
		mv $(basename $(notdir $^))_4.fastq $(basename $(notdir $^))_2.fastq
	fi


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


### Assembly Alignment:

align-assembly: $(SRA-alignments)

$(DATA_DIR)/%.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/%_megahit/final.contigs.fa
	cd $(DATA_DIR)
	nucmer --delta=$@ $^



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
