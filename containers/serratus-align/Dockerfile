FROM serratus-base:latest AS build_base


# Software Information
# ENV SAMTOOLSVERSION='1.10' # from serratus-base
ENV BOWTIEVERSION='2.4.1'
ENV GATKVERSION='4.1.5.0'
ENV BWAVERSION='0.7.17'
#ENV MBVERSION='1.5.0'


#==========================================================
# Dependencies ============================================
#==========================================================
# # Libraries for bowtie2
#RUN sudo yum -y install perl-devel libtbb-dev@testing

#==========================================================
# Install Software ========================================
#==========================================================
# BOWTIE2 ======================================= 
# /usr/local/bin/bowtie2
# /usr/local/bin/bowtie2-*

# Pre-compiled binary
RUN wget --quiet https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/"$BOWTIEVERSION"/bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip &&\
  unzip bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip &&\
  rm    bowtie2-"$BOWTIEVERSION"-linux-x86_64.zip

RUN mv bowtie2-"$BOWTIEVERSION"*/bowtie2* /usr/local/bin/ &&\
  rm -rf bowtie2-"$BOWTIEVERSION"*

# Magic-BLAST ===================================
#RUN wget --quiet ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/"$MBVERSION"/ncbi-magicblast-"$MBVERSION"-x64-linux.tar.gz &&\
#  tar -xvf ncbi-magicblast-"$MBVERSION"-x64-linux.tar.gz &&\
#  rm ncbi-magicblast-"$MBVERSION"-x64-linux.tar.gz
#
# magicblast; makeblatsdb
#RUN mv ncbi-magicblast-"$MBVERSION"/bin/* /usr/local/bin/ &&\
#  rm -rf ncbi-magicblast-*

# BWA ===========================================
RUN wget https://github.com/lh3/bwa/releases/download/v"$BWAVERSION"/bwa-"$BWAVERSION".tar.bz2 &&\
  tar -xvf bwa-* && cd bwa-"$BWAVERSION" &&\
  make && mv bwa /usr/local/bin/ &&\
  rm -rf bwa*

# SeqKit ========================================
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.12.0/seqkit_linux_amd64.tar.gz &&\
  tar -xvf seqkit* && mv seqkit /usr/local/bin/ &&\
  rm seqkit_linux*

# Set up only the bowtie2 binaries we need
RUN mkdir -p /opt/bowtie2-align \
 && cp /usr/local/bin/{bowtie2,bowtie2-align-s} /opt/bowtie2-align

FROM amazonlinux:2 AS runtime

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='aligner'
ARG VERSION='0.1.4'

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, bowtie2, magicblast, makeblatsdb, seqkit"

RUN yum -y install python3 perl \
 && alias python=python3 \
 && curl -O https://bootstrap.pypa.io/get-pip.py \
 && python3 get-pip.py \
 && rm get-pip.py \
 && pip install boto3 awscli \
 && yum clean all

# bowtie2
COPY --from=build_base /opt/bowtie2-align/* /usr/local/bin/

# samtools
COPY --from=build_base /usr/local/bin/samtools /usr/local/bin/

#==========================================================
# Serratus Initialize =====================================
#==========================================================
WORKDIR /home/serratus
COPY worker.sh ./
COPY serratus-align/*sh ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "align", "./serratus-align.sh"]
