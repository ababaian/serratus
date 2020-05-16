### Automating the comparison of Kraken, kallisto, and snap-express

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
DATA_DIR ?= /media/storage/assembly
OUT_DIR  ?= $(CURDIR)/
INSTALL_DIR ?= $(CURDIR)/third-party

LD_LIBRARY_PATH = /usr/local/lib
export LD_LIBRARY_PATH

PATH := $(PATH):$(INSTALL_DIR)/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$(INSTALL_DIR)/ncbi-blast-2.10.0+/bin
export PATH


#### Targets:

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


### Stage Data:
$(DATA_DIR)/ERR2756788_1.fastq:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	fastq-dump -I --split-files ERR2756788
	fastq-dump -I --split-files ERR2756787
	fastq-dump -I --split-files ERR3569452


$(DATA_DIR)/cov2m.fa:
	aws s3 cp s3://serratus-public/seq/cov2m/cov2m.fa

$(HOME)/.ncbirc:
	@echo "[BLAST]" > .ncbirc
	@echo "BLASTDB=$(DATA_DIR)/blastdbs" >> .ncbirc

update-blastdbs:
	mkdir -p $(DATA_DIR)/../blastdbs
	cd $(DATA_DIR)/../blastdbs
	date; time update_blastdb.pl --decompress nt; echo $$?; date

### Run Assembly:
## ~30 minutes on 72 cores, 144G RAM, and 200G instance swap
run-assembly: $(DATA_DIR)/ERR2756787_megahit/final.contigs.fa \
		$(DATA_DIR)/ERR3569452_megahit/final.contigs.fa \
		$(DATA_DIR)/ERR2756788_megahit/final.contigs.fa

$(DATA_DIR)/ERR2756788_megahit/final.contigs.fa: $(DATA_DIR)/ERR2756788_1.fastq
	megahit -1 $(DATA_DIR)/ERR2756788_1.fastq \
		-2 $(DATA_DIR)/ERR2756788_2.fastq \
		-o $(DATA_DIR)/megahit_out \
		--mem-flag 2

$(DATA_DIR)/ERR2756787_megahit/final.contigs.fa: $(DATA_DIR)/ERR2756787_1.fastq
	megahit -1 $(DATA_DIR)/ERR2756787_1.fastq \
		-2 $(DATA_DIR)/ERR2756787_2.fastq \
		-o $(DATA_DIR)/ERR2756787_megahit \
		--mem-flag 2

$(DATA_DIR)/ERR3569452_megahit/final.contigs.fa: $(DATA_DIR)/ERR3569452_1.fastq
	megahit -1 $(DATA_DIR)/ERR3569452_1.fastq \
		-2 $(DATA_DIR)/ERR3569452_2.fastq \
		-o $(DATA_DIR)/ERR3569452_megahit \
		--mem-flag 2

align-assembly: $(DATA_DIR)/ERR3569452.delta \
	$(DATA_DIR)/ERR2756787.delta \
	$(DATA_DIR)/ERR2756788.delta

$(DATA_DIR)/ERR3569452.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR3569452_megahit/final.contigs.fa
	cd $(DATA_DIR)
	nucmer --prefix=ERR3569452 $^

$(DATA_DIR)/ERR2756787.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR2756787_megahit/final.contigs.fa
	cd $(DATA_DIR)
	nucmer --prefix=ERR2756787 $^

$(DATA_DIR)/ERR2756788.delta: $(DATA_DIR)/cov2m.fa $(DATA_DIR)/ERR2756788_megahit/final.contigs.fa
	cd $(DATA_DIR)
	nucmer --prefix=ERR2756788 $^

## This shows that k141_146741 wholly contains three fragments in cov2m, and is ~29k:
show-coords:
	cd $(DATA_DIR)
	show-coords -c -q -g -o -l -T ERR2756788.delta > ERR2756788-alignments.txt
	show-coords -c -q -g -o -l -T ERR2756787.delta > ERR2756787-alignments.txt
	show-coords -c -q -g -o -l -T ERR3569452.delta > ERR3569452-alignments.txt

stage-results-to-s3:
	cd $(DATA_DIR)
	for acc in ERR2756787 ERR2756788 ERR3569452
	do
		tar czf $${acc}_megahit.tar.gz $${acc}_megahit
		aws s3 cp $${acc}_megahit.tar.gz s3://serratus-public/notebook/200515_ta/
	done
