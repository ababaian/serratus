#!/bin/sh
#
# makeDownloader v0.1.2
# Script to make "serratus-Downloader" AMI
#
# Base Image: Amazon Linux 2
# AMI: ami-0e8c04af2729ff1bb
# login: ec2-user@<ipv4>
# base: 9 Gb
#
# Image: serratus-Downloader v0.1.2
# Desc : bioinformatics seq-database (SRA GDC) access
# AMI  : ami-0614da38b644ca19e (us-west-2)

# Software
# AWSCLI -- pre-installed on amazon linux
SAMTOOLSVERSION='1.10'
SRATOOLKITVERSION='2.10.4'
GDCVERSION='1.5.0'
PICARDVERSION='2.22.0'
GATKVERSION='4.1.5.0'
#LIBMAUSVERSION='2.0.705'
#BIOBAMBAMVERSION='2.0.157'

# DEPENDENCY ====================================
# Update core
sudo yum update
sudo yum clean all
sudo yum install git

# Python3 3.7.4 and pip3
sudo yum install python3
sudo yum install python3-devel
alias python=python3

curl -O https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user
rm get-pip.py

# Libraries for htslib
sudo yum install make gcc
sudo yum install unzip bzip2 bzip2-devel xz-devel zlib-devel
sudo yum install curl-devel
sudo yum install ncurses-devel openssl-devel

# Libraries/software for GATK/JAVA
sudo yum install java

# Assorted libs for scripts
sudo yum install pigz

# Libraries for biobambam2
#sudo yum install gcc-c++

# SAMTOOLS ======================================
# /usr/local/bin/samtools
wget -O samtools-"$SAMTOOLSVERSION".tar.bz2 \
  https://github.com/samtools/samtools/releases/download/"$SAMTOOLSVERSION"/samtools-"$SAMTOOLSVERSION".tar.bz2

tar xvjf samtools-"$SAMTOOLSVERSION".tar.bz2 && rm samtools-"$SAMTOOLSVERSION".tar.bz2

cd samtools-"$SAMTOOLSVERSION"
  make
  sudo make install
cd .. && rm -rf samtools-*

# SRATOOLKIT=====================================
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/"$SRATOOLKITVERSION"/setup-apt.sh
sudo ./setup-apt.sh && rm ./setup-apt.sh
# Add to $PATH:/usr/local/ncbi/sra-tools/bin
source /etc/profile.d/sra-tools.sh

# Test command:
# fastq-dump --stdout -X 2 SRR390728

# See: https://github.com/ncbi/sra-tools/wiki/04.-Cloud-Credentials
# To access SRA cloud-data, you'll need to provide 
# your AWS (Amazon Web Services) access key or GCP 
# (Google Cloud Platform) service account to vdb-config.

# GDC-CLIENT ====================================
wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v"$GDCVERSION"_Ubuntu_x64.zip

unzip gdc-client_v"$GDCVERSION"_Ubuntu_x64.zip
rm    gdc-client_v"$GDCVERSION"_Ubuntu_x64.zip
sudo mv gdc-client /usr/local/bin/

# Clean-up
sudo yum clean all
sudo rm -rf /var/cache/yum

# Save AMI
# ami (us-west-2): ami-059b454759561d9f4

# PICARD ========================================
wget https://github.com/broadinstitute/picard/releases/download/"$PICARDVERSION"/picard.jar
chmod 755 picard.jar

# ENV VARIABLE
PICARD='/usr/local/bin/picard.jar'
sudo mv picard.jar $PICARD

# Default: picard -h
alias picard="java -jar $PICARD"
# Java optimized
# "java -jvm-args -jar $PICARD -h"

# GATK ==========================================
wget https://github.com/broadinstitute/gatk/releases/download/"$GATKVERSION"/gatk-"$GATKVERSION".zip
unzip gatk*
rm gatk*zip

sudo mv gatk-$GATKVERSION /usr/local/share/gatk
sudo ln -s /usr/local/share/gatk/gatk /usr/local/bin/

# # BIOBAMBAM2 ====================================
# # DEPENDENCY HELL -- ABANDON HOPE ALL YE WHO ENTER HERE
# # libmaus2.0.705 dependency
# wget https://gitlab.com/german.tischler/libmaus2/-/archive/2.0.705-release-20200303192236/libmaus2-2.0.705-release-20200303192236.zip
# unzip libmaus2*

# cd libmaus2*
#   sudo mkdir /usr/local/libmaus2
#   sudo chown ec2-user /usr/local/libmaus2

#   ./configure --prefix=/usr/local/libmaus2
#   make 
#   make install
# cd ..

# # biobambam2
# wget https://gitlab.com/german.tischler/biobambam2/-/archive/2.0.156-release-20200224081934/biobambam2-2.0.156-release-20200224081934.zip
# unzip biobambam2*

# cd biobambam2*
#   sudo mkdir /usr/local/biobambam2
#   sudo chown ec2-user /usr/local/biobambam2

#   sudo ./configure --with-libmaus2=/usr/local/libmaus2/2.0.705
# make install

