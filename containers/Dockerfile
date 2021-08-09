#==========================================================
# SERRATUS-BASE CONTAINER =================================
#==========================================================
# Docker Base: amazon linux2
FROM amazonlinux:2 AS serratus-base

## Build/test container for serratus
# sudo yum install -y docker git
# sudo service docker start
# git clone https://github.com/ababaian/serratus.git; cd serratus
# sudo docker build -t serratus-base:0 -t serratus-base:latest -f docker/Dockerfile .
# sudo docker build -t serratus-align:0 -t serratus-align:latest -f docker/serratus-align/Dockerfile .
# sudo docker run --rm --entrypoint /bin/bash -it serratus-align:0

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG PROJECT='serratus'
ARG TYPE='base'
ARG VERSION='0.3.0'

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
RUN yum -y install tar wget gzip which sudo shadow-utils \
           util-linux byacc

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

# Increase the default chunksize for `aws s3 cp`.  By default it is 8MB,
# which results in a very high number of PUT and POST requests.  These
# numbers have NOT been experimented on, but chosen to be just below the
# max size for a single-part upload (5GB).  I haven't pushed it higher
# because I don't want to test the edge cases where a filesize is around
# the part limit.
RUN aws configure set default.s3.multipart_threshold 4GB \
 && aws configure set default.s3.multipart_chunksize 4GB

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
# N/A
