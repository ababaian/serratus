FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='aligner'
ARG VERSION='0.1.0'

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
# # Libraries for bowtie2
#RUN sudo yum -y install perl-devel libtbb-dev@testing

#==========================================================
# Install Software ========================================
#==========================================================
# BOWTIE2 ======================================= 
# /usr/local/bin/bowtie2
# /usr/local/bin/bowtie2-*

# install git
RUN yum -y install git

#compile bowtie from sources
RUN git clone --recursive https://github.com/BenLangmead/bowtie2.git &&\
  cd bowtie2 &&\
  make POPCNT_CAPABILITY=0 NO_TBB=1 &&\
  mv bowtie2 /usr/local/bin/ &&\
  mv bowtie2-* /usr/local/bin/ &&\
  cd .. &&\
  rm -rf bowtie2

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
