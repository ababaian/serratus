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
DATA_DIR ?= /media/storage/kraken2
OUT_DIR  ?= $(CURDIR)/

PATH := $(PATH):/usr/local/ncbi/sra-tools/bin
export PATH

#### Targets:

### Install Dependencies:
install-packages:
	sudo apt-get install libxml-libxml-perl curl

install-sra-toolkit:
	mkdir -p ../third-party
	cd ../third-party
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.5/setup-apt.sh
	chmod 755 ./setup-apt.sh
	sudo ./setup-apt.sh


### Data:
stage-samples:
	mkdir -p $(DATA_DIR)/samples
	aws s3 cp s3://sra-pub-run-2/SRR11728619/SRR11728619.2 $(DATA_DIR)/samples/

### Run Kraken2 on ERR2756788 'Frankie':
run-kraken2-frankie:
	date; time ~/repos/snap-express/third-party/bin/kraken2 --paired --report profile.kraken2 --use-names --db /media/storage/kraken2-db --threads `nproc` /media/storage/assembly/ERR2756788_1.fastq /media/storage/assembly/ERR2756788_2.fastq > classified-reads.kraken2.txt ; echo $?; date
