# Docker Base: amazon linux2
FROM amazonlinux:2 AS build_base

# docker build --no-cache -t serratus-align:0 .
# docker run --rm --entrypoint /bin/bash -it serratus-align:0

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG PROJECT='serratus'
ARG TYPE='aligner'
ARG VERSION='0.1'

# Software Information (pass to shell)
ENV SAMTOOLSVERSION='1.10'
ENV BOWTIEVERSION='2.4.1'
ENV GATKVERSION='4.1.5.0'

# Additional Metadata
LABEL container.base.image="amazonlinux:2"
LABEL project.name=${SOFTWARE}
LABEL project.website="https://github.com/ababaian/serratus"
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, bowtie2, gatk"

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

# # Libraries for bowtie2
#RUN sudo yum -y install perl-devel libtbb-dev@testing


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

# BOWTIE2 ======================================= 
# /usr/local/bin/bowtie2
# /usr/local/bin/bowtie2-*

# Pre-compiled binary
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/"$BOWTIEVERSION"/bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip &&\
  unzip bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip &&\
  rm    bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip

RUN mv bowtie2-"$BOWTIEVERSION"*/bowtie2* /usr/local/bin/ 
RUN rm -rf bowtie2-"$BOWTIEVERSION"*

#==========================================================
# Serratus Initialize =====================================
#==========================================================
# Set working directory to /home/serratus
# Create a serratus user fro use of FUSE / fusera
RUN adduser serratus
ARG BASEDIR=/home/serratus
ENV BASEDIR=/home/serratus
RUN echo 'alias serratus="sudo -u serratus"' >> /root/.bashrc

RUN yum install -y util-linux # flock

WORKDIR /home/serratus
COPY scripts/worker.sh ./
COPY scripts/serratus-align/* ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
#Usage: docker exec serratus-align -u SCHED [OPTIONS]

ENTRYPOINT ["./worker.sh", "align", "./serratus-align.sh"]

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
