# Docker Base: amazon linux2
FROM amazonlinux:2 AS build_base

# Deploy / Development
# AMI: Amazon Linux 2 with extras
# EC2 AMI: ami-079f731edfe27c29c
# sudo yum -y install docker
# sudo service docker start

# sudo docker build --no-cache -t serratus-dl:0 -t serratus-dl:latest .
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
ARG TYPE='downloader'
ARG VERSION='0.1.2'

# Software Information (pass to shell)
ENV SAMTOOLSVERSION='1.10'
ENV SRATOOLKITVERSION='2.10.4'
ENV GDCVERSION='1.5.0'
ENV PICARDVERSION='2.22.0'
#ARG GATKVERSION='4.1.5.0'

# Additional Metadata
LABEL container.base.image="amazonlinux:2"
LABEL project.name=${SOFTWARE}
LABEL project.website="https://github.com/ababaian/serratus"
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: downloader and fq-splitter"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, sra-toolkit, gdc-client"

#==========================================================
# Dependencies ============================================
#==========================================================
# For development only
RUN yum -y install vim htop

# Update Core
RUN yum -y update
RUN yum -y install tar wget gzip which

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

# Libraries / User for FUSE
RUN yum -y install sudo shadow-utils \
    fuse fuse-devel

# Libraries for serratus-downloader scripts
RUN yum -y install pigz

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

# SRATOOLKIT ====================================
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKITVERSION}/sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    tar xzf sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    rm -f sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    mv ./sratoolkit.${SRATOOLKITVERSION}-centos_linux64 /opt &&\
    echo -e "export PATH=/opt/sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin:\$PATH" >\
            /etc/profile.d/sratoolkit.sh &&\
    chmod 755 /etc/profile.d/sratoolkit.sh &&\
    mkdir /etc/ncbi

ENV PATH "$PATH:/opt/sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin"

# GDC-CLIENT ===================================
 RUN wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    unzip gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    rm    gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    mv gdc-client /usr/local/bin/

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
COPY scripts/serratus-dl/* ./
COPY docker/serratus-dl/VDB_user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN mkdir -p /root/.ncbi; vdb-config --report-cloud-identity yes

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
#Usage: docker run serratus-dl -u $SCHEDULER:$PORT [OPTIONS]

ENTRYPOINT ["./serratus-dl.sh"]

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
