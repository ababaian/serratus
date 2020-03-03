#!/bin/sh
#
# makeAlinger v0.1
# Script to make "serratus-Aligner" AMI
#
# Base Image: Amazon Linux 2
# AMI: ami-0e8c04af2729ff1bb
# login: ec2-user@<ipv4>
# base: 9 Gb
#
# Image: serratus-Aligner
# Desc : (v0.1) bioinformatics alignment - bowtie2
# AMI  : ami-059b454759561d9f4 (us-west-2)

# Software
SAMTOOLSVERSION='1.10'
BOWTIEVERSION='2.4.1'
#GATK=''

# DEPENDENCY ====================================
# Update core
sudo yum update
sudo yum clean all

# Python3 3.7.4 and pip3
sudo yum install python3
sudo yum install python3-devel
alias python=python3

curl -O https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user
rm get-pip.py

# Libraries for htslib
sudo yum install make gcc libc-dev
sudo yum install unzip bzip2-devel xz-devel zlib-devel
sudo yum install ncurses-devel
sudo yum install curl-devel

# # Libraries for bowtie2
# sudo apk add g++
# sudo apk add perl perl-dev
# sudo apk add libtbb-dev@testing

# SAMTOOLS ======================================
# /usr/local/bin/samtools
wget -O samtools-"$SAMTOOLSVERSION".tar.bz2 \
  https://github.com/samtools/samtools/releases/download/"$SAMTOOLSVERSION"/samtools-"$SAMTOOLSVERSION".tar.bz2

tar xvjf samtools-"$SAMTOOLSVERSION".tar.bz2 && rm samtools-"$SAMTOOLSVERSION".tar.bz2

cd samtools-"$SAMTOOLSVERSION"
  make
  sudo make install
cd .. && rm -rf samtools-*

# bowtie2 =======================================
# /usr/local/bin/bowtie2
# /usr/local/bin/bowtie2-*

# Pre-compiled binary
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/"$BOWTIEVERSION"/bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip
unzip bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip
rm    bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip

sudo mv bowtie2-*/bowtie2* /usr/local/bin/
rm bowtie2*

# Save AMI