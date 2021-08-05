FROM serratus-base:latest AS build_base

# Software Information
# ENV SAMTOOLSVERSION='1.10' # from serratus-base
ENV BOWTIEVERSION='2.4.1'
ENV GATKVERSION='4.1.5.0'
ENV BWAVERSION='0.7.17'
ENV DIAMONDVERSION='2.0.6-dev'
ENV INFERNALVERSION='1.1.4'
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

# DIAMOND ======================================
# RUN wget --quiet https://github.com/bbuchfink/diamond/releases/download/v"$DIAMONDVERSION"/diamond-linux64.tar.gz &&\
#   tar -xvf diamond-linux64.tar.gz &&\
#   rm    diamond-linux64.tar.gz &&\
#   mv    diamond /usr/local/bin/

# Use serratus-built dev version
RUN wget --quiet https://serratus-public.s3.amazonaws.com/bin/diamond &&\
    chmod 755 diamond &&\
    mv    diamond /usr/local/bin/

# BWA ===========================================
RUN wget https://github.com/lh3/bwa/releases/download/v"$BWAVERSION"/bwa-"$BWAVERSION".tar.bz2 &&\
  tar -xvf bwa-* && cd bwa-"$BWAVERSION" &&\
  make && mv bwa /usr/local/bin/ &&\
  rm -rf bwa*

# SeqKit ========================================
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.12.0/seqkit_linux_amd64.tar.gz &&\
  tar -xvf seqkit* && mv seqkit /usr/local/bin/ &&\
  rm seqkit_linux*

# Bowtie2 =======================================
# Set up only the bowtie2 binaries we need
RUN mkdir -p /opt/bowtie2-align \
 && cp /usr/local/bin/{bowtie2,bowtie2-align-s} /opt/bowtie2-align

# INFERNAL ======================================
# cmsearch is required alone (in theory)
# build from source
RUN wget http://eddylab.org/infernal/infernal-"$INFERNALVERSION"-linux-intel-gcc.tar.gz &&\
  tar -xvf infernal-1.1.4-linux-intel-gcc.tar.gz &&\
  cd infernal* &&\
  bash configure && make && make install &&\
  cd .. && rm -rf infernal*

# binaries copy (does not work)
#  chmod 755 infernal*/binaries/* &&\
#  cp infernal*/binaries/* /usr/local/bin/ &&\
#  rm -rf infernal*

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='aligner'
ARG VERSION='0.3.0'

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, bowtie2, magicblast, makeblatsdb, seqkit"

# Environment setup
ENV BASEDIR=/home/serratus

RUN yum -y install perl jq unzip \
 && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
 && unzip awscliv2.zip \
 && ./aws/install \
 && rm awscliv2.zip \
 && rm -r ./aws \
 && yum -y remove unzip \
 && yum clean all \
#   aws configuration
 && aws configure set default.s3.multipart_threshold 4GB \
 && aws configure set default.s3.multipart_chunksize 4GB

# bowtie2
COPY --from=build_base /opt/bowtie2-align/* /usr/local/bin/

# diamond
COPY --from=build_base /usr/local/bin/diamond /usr/local/bin/

# samtools
COPY --from=build_base /usr/local/bin/samtools /usr/local/bin/

# INFERNAL
COPY --from=build_base /usr/local/bin/cmsearch /usr/local/bin/
#COPY --from=build_base /usr/local/bin/cmscan /usr/local/bin/
# ...

# flock
COPY --from=build_base /usr/bin/flock /usr/bin/


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
