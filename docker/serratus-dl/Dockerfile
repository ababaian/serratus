FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='downloader'
ARG VERSION='0.1.0'

# Software Information
# ENV SAMTOOLSVERSION='1.10'  # from serratus-base
ENV SRATOOLKITVERSION='2.10.4'
ENV GDCVERSION='1.5.0'
#ENV PICARDVERSION='2.22.0'

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: downloader and fq-splitter"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, sra-toolkit, gdc-client"

#==========================================================
# Dependencies ============================================
#==========================================================
# Libraries for serratus-downloader scripts
RUN yum -y install pigz

#==========================================================
# Install Software ========================================
#==========================================================
# SRATOOLKIT ====================================
RUN wget --quiet https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKITVERSION}/sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    tar xzf sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    rm -f sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz &&\
    mv ./sratoolkit.${SRATOOLKITVERSION}-centos_linux64 /opt &&\
    echo -e "export PATH=/opt/sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin:\$PATH" >\
            /etc/profile.d/sratoolkit.sh &&\
    chmod 755 /etc/profile.d/sratoolkit.sh &&\
    mkdir /etc/ncbi

ENV PATH "$PATH:/opt/sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin"

# GDC-CLIENT ===================================
 RUN wget --quiet https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    unzip gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    rm    gdc-client_v${GDCVERSION}_Ubuntu_x64.zip &&\
    mv gdc-client /usr/local/bin/

#==========================================================
# Serratus Initialize =====================================
#==========================================================
# Set working directory to /home/serratus
# Create a serratus user fro use of FUSE / fusera

WORKDIR /home/serratus
COPY scripts/worker.sh ./
COPY scripts/serratus-dl/* ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "split", "./serratus-dl.sh"]
