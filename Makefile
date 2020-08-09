#### Assumptions:
## * AWS CLI installed & configured
## * sudo apt-get install make python2.7

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


### Local definitions:

export PATH := $(PATH):$(CURDIR)/third-party/mash-Linux64-v2.2:$(CURDIR)/third-party/sratoolkit.2.10.5-ubuntu64/bin


#### Targets:

install-mash:
	mkdir -p third-party
	cd third-party
	wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar
	tar xf mash-Linux64-v2.2.tar
	cd mash-Linux64-v2.2.tar

install-sra:
	mkdir -p third-party
	cd third-party
	wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
	tar xzf sratoolkit.tar.gz

download-mash-test-data:
	mkdir -p test
	cd test
	fasterq-dump --split-files SRR11454614
	cat SRR11180057*.fastq > SRR11454614.fastq
	fasterq-dump --split-files SRR9658384
	cat SRR9658384*.fastq > SRR9658384.fastq
	aws s3 cp s3://serratus-public/seq/cov2r/cov2.fa.gz .
	aws s3 sync s3://serratus-public/notebook/200411/fq fq

test-mash:
	cd test
	@echo "Sketching the cov2 DB":
	time mash sketch cov2.fa.gz 
	@echo "Testing cov2 vs. SRR9658384 reads (should be zero):"
	time mash screen -p 2 cov2.fa.gz.msh SRR9658384.fastq
	@echo "Testing cov2 vs. SRR11454614 reads (should be large):"
	time mash screen -p 2 cov2.fa.gz.msh SRR11454614.fastq


test-mash-sensitivity:
	cd test
	mash sketch -k 15 -s 5000 cov2.fa.gz
	for sim in 0 30 300 1500 3000 4500 6000 7500 9000 10500 12000
	do
		echo Simulation $$sim
		zcat fq/sim.cov.$${sim}_1.fq fq/sim.cov.$${sim}_2.fq.gz \
			> sim.cov.$${sim}.fq
		time mash screen -p 2 cov2.fa.gz.msh sim.cov.$${sim}.fq 
	done


### Conda

## I install it to /home/ubuntu/.miniconda3
install-miniconda:
	wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	if [ `sha256sum Miniconda3-latest-Linux-x86_64.sh | cut -d' ' -f 1` == 957d2f0f0701c3d1335e3b39f235d197837ad69a944fa6f5d8ad2c686b69df3b ]
	then
		chmod 700 Miniconda3-latest-Linux-x86_64.sh
		./Miniconda3-latest-Linux-x86_64.sh
	fi

#
# To activate this environment, use
#
#     $ conda activate snakemake
#
# To deactivate an active environment, use
#
#     $ conda deactivate

install-snakemake:
	conda install -c conda-forge mamba
	mamba create -c conda-forge -c bioconda -n snakemake snakemake

install-igv:
	cd third-party
	wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_Linux_2.8.2.zip


### Utilities:
## Use lsblk to find a suitable ephemeral drive to use, then call it like:
## SWAP_DEVICE=/dev/nvme1n1 SWAP_SIZE=400G make set-up-swap
storage-dir := /media/storage
set-up-swap:
	sudo mkdir -p $(storage-dir)
	sudo mkfs.ext4 $(SWAP_DEVICE)
	sudo mount $(SWAP_DEVICE) $(storage-dir)
	sudo fallocate -l $(SWAP_SIZE) $(storage-dir)/swapfile
	sudo chmod 600 $(storage-dir)/swapfile
	sudo mkswap $(storage-dir)/swapfile
	sudo swapon $(storage-dir)/swapfile
