# Docker Base: amazon linux2
FROM amazonlinux:2 AS build_base

# docker build --no-cache -t serratus-merge:0 .
# docker run --rm --entrypoint /bin/bash -it serratus-merge:0

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG PROJECT='serratus'
ARG TYPE='merge'
ARG VERSION='0.1'

# Software Information (pass to shell)
ENV SAMTOOLSVERSION='1.10'

# Additional Metadata
LABEL container.base.image="amazonlinux:2"
LABEL project.name=${SOFTWARE}
LABEL project.website="https://github.com/ababaian/serratus"
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools"

#==========================================================
# Dependencies ============================================
#==========================================================
# For development only
RUN yum -y install vim htop

# Update Core
RUN yum -y update
RUN yum -y install tar wget gzip which

# User libraries
RUN yum -y install sudo shadow-utils

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
RUN cd $BASEDIR; mkdir -p scripts

# Serratus Scripts ==============================
# URL to latest serratus repo (raw)
ARG GIT_SERRATUS="https://raw.githubusercontent.com/ababaian/serratus/master/"

## /home/serratus/serratus-merge.sh
RUN wget -O $BASEDIR/serratus-merge.sh \
  $GIT_SERRATUS/scripts/serratus-merge/serratus-merge.sh

## /home/serratus/scripts/
  # run_download.sh
RUN wget -O $BASEDIR/scripts/run_merge.sh \
  $GIT_SERRATUS/scripts/serratus-merge/run_merge.sh

# Set script permissions
RUN chmod 755 $BASEDIR/serratus-merge.sh
RUN chmod 755 $BASEDIR/scripts/*

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
#Usage: docker exec serratus-align -u SCHED [OPTIONS]

ENTRYPOINT ["/bin/bash", "/home/serratus/serratus-merge.sh"]

#==========================================================
# DEV  ====================================================
#==========================================================
## EC2
# sudo yum -y install docker
# sudo service docker start

# sudo docker build -t serratus-dl:0 .
# sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0

# sudo docker run --rm serratus-dl:0
# sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0