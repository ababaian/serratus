FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='aligner'
ARG VERSION='0.1.2'

# Software Information
# ENV SAMTOOLSVERSION='1.10' # from serratus-base
ENV BOWTIEVERSION='2.4.1'
ENV GATKVERSION='4.1.5.0'

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, bowtie2, gatk"

#==========================================================
# Dependencies ============================================
#==========================================================
RUN yum install -y util-linux # flock

# # Libraries for bowtie2
#RUN sudo yum -y install perl-devel libtbb-dev@testing

#==========================================================
# Install Software ========================================
#==========================================================
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
WORKDIR /home/serratus
COPY scripts/worker.sh ./
COPY scripts/serratus-align/* ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "align", "./serratus-align.sh"]
