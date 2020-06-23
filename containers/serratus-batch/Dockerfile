FROM serratus-base:latest AS build_base

# Versions
ENV SRATOOLKITVERSION='2.10.4'
ENV BOWTIEVERSION='2.4.1'

# Bowtie2 - download and install
RUN wget -O bowtie2.zip --quiet "https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIEVERSION}/bowtie2-${BOWTIEVERSION}-linux-x86_64.zip" \
 && unzip bowtie2.zip \
 && rm bowtie2.zip \
 && mkdir /opt/bowtie2-align \
 && mv bowtie2*/{bowtie2,bowtie2-align-s} /opt/bowtie2-align \
 && rm -rf bowtie2*

## SRAToolkit
RUN wget --quiet "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKITVERSION}/sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz" \
 && tar xzf "sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz" \
 && rm -f "sratoolkit.${SRATOOLKITVERSION}-centos_linux64.tar.gz" \
 && mkdir -p /opt/sratools \
    # Keep sratools grouped together, so its easy to copy them all out into the runtime
 && mv ./sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin/{vdb-config*,prefetch*,fastq-dump*,fasterq-dump*,sratools*} /opt/sratools \
    # Install into /usr/local/bin for the rest of the build
 && cp -r /opt/sratools/* /usr/local/bin \
 && mkdir /etc/ncbi

# Copy in config for sra tools
COPY serratus-dl/VDB_user-settings.mkfg /root/.ncbi/user-settings.mkfg
RUN vdb-config --report-cloud-identity yes

FROM amazonlinux:2 AS runtime

LABEL container.description="serratus: end-to-end batch container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools, bowtie2, sratoolkit"

# aws cli, plus dependencies
# -merge has its own python dependency, so do a full python/pip install
RUN yum -y install python3 perl \
 && alias python=python3 \
 && curl -O https://bootstrap.pypa.io/get-pip.py \
 && python3 get-pip.py \
 && rm get-pip.py \
 && pip install boto3 awscli \
 && yum clean all \
#   aws configuration
 && aws configure set default.s3.multipart_threshold 4GB \
 && aws configure set default.s3.multipart_chunksize 4GB

# bowtie2
COPY --from=build_base /opt/bowtie2-align/* /usr/local/bin/

# samtools
COPY --from=build_base /usr/local/bin/samtools /usr/local/bin/

# sratools
COPY --from=build_base /root/.ncbi /root/.ncbi
COPY --from=build_base /opt/sratools/ /usr/local/bin/

# run script
WORKDIR /home/serratus
COPY ./serratus-batch/run /usr/local/bin/
COPY ./serratus-merge/serratus_summarizer.py ./summarizer.py

CMD ["run"]
