#==========================================================
# SERRATUS-BASE CONTAINER =================================
#==========================================================
# Docker Base: amazon linux2
FROM amazonlinux:2 AS serratus-base

# Development Notes:
# serratus-base
# |    |--serratus-dl
# |    |--serratus-align
# |    |     |--serratus-align-hgr1
# |    |     |--...
# |    |--serratus-merge

## Build all containers for serratus
# git clone https://github.com/ababaian/serratus.git; cd serratus
# sudo docker build -t serratus-base:0 -t serratus-base:latest -f docker/Dockerfile .
# sudo docker build -t serratus-dl:0 -t serratus-dl:latest -f docker/serratus-dl/Dockerfile .
# sudo docker build -t serratus-align:0 -t serratus-align:latest -f docker/serratus-align/Dockerfile .
# sudo docker build -t serratus-merge:0 -t serratus-merge:latest -f docker/serratus-merge/Dockerfile .

## Uploading serratus container images to AWS ECR
## Paste resulting command in terminal to authenticate
# aws ecr get-login --region us-east-1 --no-include-email

## Tag images with ECR information
# sudo docker tag serratus-base:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker tag serratus-dl:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker tag serratus-align:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker tag serratus-merge:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge

## Push images to ECR (rather host on dockerhub)
## Make sure to initialize each repository on your AWS ECR account
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge

## Run interactive serratus-dl
# sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0

# cd serratus/scheduler
# sudo docker build -t scheduler:0 .
# docker run -d --rm -p8000:8000 --name sch scheduler:0
# curl -T /path/to/SraRunInfo.csv localhost:8000/add_sra_run_info

# sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0


# Test containrs
# curl -X POST -T SraRunInfo_test.csv localhost:8000/jobs (edited) 
# curl -X POST localhost:8000/jobs/split/

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG PROJECT='serratus'
ARG TYPE='base'
ARG VERSION='0.1.2'

# Software Versions (pass to shell)
ENV SAMTOOLSVERSION='1.10'

# Additional Metadata
LABEL container.base.image="amazonlinux:2"
LABEL project.name=${PROJECT}
LABEL project.website="https://github.com/ababaian/serratus"
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus-base image"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools"

#==========================================================
# Dependencies ============================================
#==========================================================
# For development only
RUN yum -y install vim htop

# Update Core
RUN yum -y update
RUN yum -y install tar wget gzip which sudo shadow-utils

# Python3
RUN yum -y install python3 python3-devel
RUN alias python=python3
RUN curl -O https://bootstrap.pypa.io/get-pip.py &&\
    python3 get-pip.py &&\
    rm get-pip.py

# AWS S3
RUN pip install boto3 awscli
RUN yum -y install jq

# Libraries for htslib
RUN yum -y install gcc make \
    unzip bzip2 bzip2-devel xz-devel zlib-devel \
    curl-devel openssl-devel \
    ncurses-devel

#==========================================================
# Install Software ========================================
#==========================================================

# SAMTOOLS ======================================
## Download
RUN wget -O /samtools-${SAMTOOLSVERSION}.tar.bz2 \
  https://github.com/samtools/samtools/releases/download/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 &&\
  tar xvjf /samtools-${SAMTOOLSVERSION}.tar.bz2 && rm /samtools-${SAMTOOLSVERSION}.tar.bz2

# Install
RUN cd samtools-${SAMTOOLSVERSION} && make && make install

# Cleanup
RUN cd .. && rm -rf samtools-${SAMTOOLSVERSION}

#==========================================================
# Serratus Initialize =====================================
#==========================================================
# Set working directory to /home/serratus
# Create a serratus user fro use of FUSE / fusera
RUN adduser serratus
ARG BASEDIR=/home/serratus
ENV BASEDIR=/home/serratus
RUN echo 'alias serratus="sudo -u serratus"' >> /root/.bashrc

WORKDIR /home/serratus

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
# N/A