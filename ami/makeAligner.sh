#!/bin/sh
#
# makeAlinger
# Script to make "Aligner" AMI
#
# Base Image: Alpine Linux 3.11
# AMI: ami-050dd0423825ae4cd
# login: alpine@<ipv4>
#

# Software
# - samtools
SAMTOOLSVERSION='1.10'
# - bowtie2
BOWTIEVERSION='2.4.1'
# - GATK

# DEPENDENCY ====================================
# Update core
sudo apk update
sudo apk upgrade                          # 96 Mb
sudo apk add bash                         # 97 Mb

## Add edge repositories (required for libtpp)
## to /etc/apk/repositories
# http://dl-cdn.alpinelinux.org/alpine/v3.11/main
# http://dl-cdn.alpinelinux.org/alpine/v3.11/community
# @edge http://dl-cdn.alpinelinux.org/alpine/edge/main
# @edgecommunity http://dl-cdn.alpinelinux.org/alpine/edge/community
# @testing http://dl-cdn.alpinelinux.org/alpine/edge/testing

# Python3
sudo apk add python3                      #155 Mb
sudo apk add py-pip                       #208 Mb

# AWS CLI
sudo pip install boto3 awscli

# Libraries for htslib
sudo apk add make gcc libc-dev               #315 Mb
sudo apk add unzip bzip2-dev xz-dev zlib-dev #316 Mb
sudo apk add ncurses-dev                     #316 Mb
sudo apk add curl-dev                        #319 Mb

# # Libraries for bowtie2
sudo apk add g++                             #379 Mb
sudo apk add perl perl-dev                   #423 Mb
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

sudo wget -O tmp.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-linux-x86_64.zip