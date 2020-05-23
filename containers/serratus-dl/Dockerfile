FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='downloader'
ARG VERSION='0.1.4'

# Software Information
# ENV SAMTOOLSVERSION='1.10'  # from serratus-base
ENV SRATOOLKITVERSION='2.10.4'
#ENV GDCVERSION='1.5.0'
#ENV BVVERSION='0.1'
#ENV PICARDVERSION='2.22.0'

#==========================================================
# Dependencies ============================================
#==========================================================
# Libraries for serratus-downloader scripts
#RUN yum -y install pigz

# EPEL + Parallel
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
RUN yum install -y parallel procps-ng; echo will cite | parallel --citation 1>/dev/null 2>/dev/null

#==========================================================
# Install Software ========================================
#==========================================================
# SRATOOLKIT ====================================
RUN wget --quiet https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKITVERSION}/sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    tar xzf sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    rm -f sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    mkdir -p /opt/sratools &&\
    # Keep sratools grouped together, so its easy to copy them all out into the runtime
    mv ./sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin/{vdb-config*,prefetch*,fastq-dump*,fasterq-dump*,sratools*} /opt/sratools &&\
    # Install into /usr/local/bin for the rest of the build
    cp -r /opt/sratools/* /usr/local/bin &&\
    mkdir /etc/ncbi

# GDC-CLIENT ===================================
# RUN wget --quiet https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
#    unzip gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
#    rm    gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
#    mv gdc-client /usr/local/bin/

# BVFILTER =====================================
#COPY bin/bvfilter /usr/local/bin/

#==========================================================
# Serratus Initialize =====================================
#==========================================================
# Set working directory to /home/serratus
# Create a serratus user fro use of FUSE / fusera
COPY serratus-dl/VDB_user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN mkdir -p /root/.ncbi; vdb-config --report-cloud-identity yes

FROM amazonlinux:2 AS runtime

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: downloader and fq-splitter"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, sra-toolkit"

# Environment setup
ENV BASEDIR=/home/serratus

# aws cli, plus dependencies
RUN yum -y install ps jq python3 perl \
 && alias python=python3 \
 && curl -O https://bootstrap.pypa.io/get-pip.py \
 && python3 get-pip.py \
 && rm get-pip.py \
 && pip install boto3 awscli \
 && yum clean all \
#   aws configuration
 && aws configure set default.s3.multipart_threshold 4GB \
 && aws configure set default.s3.multipart_chunksize 4GB

# sratools
COPY --from=build_base /root/.ncbi /root/.ncbi
COPY --from=build_base /opt/sratools/ /usr/local/bin/

# parallel
COPY --from=build_base /usr/bin/parallel /usr/bin/
# echo will cite | parallel --citation 1>/dev/null 2>/dev/null

# flock
COPY --from=build_base /usr/bin/flock /usr/bin/

WORKDIR /home/serratus
COPY worker.sh ./
COPY serratus-dl/*sh ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "split", "./serratus-dl.sh"]
