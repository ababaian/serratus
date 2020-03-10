# Docker Base: amazon linux2
FROM amazonlinux:2 AS build_base

# Deploy / Development
# AMI: Amazon Linux 2 with extras
# EC2 AMI: ami-079f731edfe27c29c
# sudo yum -y install docker
# sudo service docker start

# sudo docker build -t serratus-dl:0 .
# sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0

# Run container with /dev/fuse access point
# sudo docker run --rm --entrypoint /bin/bash -it --privileged --cap-add SYS_ADMIN --cap-add MKNOD --device /dev/fuse serratus-dl:0

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

RUN echo 'export PATH=/opt/sratoolkit.${SRATOOLKITVERSION}-centos_linux64/bin:"$PATH"' >> /root/.bashrc

# FUSERA ========================================
# # installed as user `serratus`
# # fusera cannot be run as root
# RUN wget https://raw.githubusercontent.com/mitre/fusera/master/install.sh &&\
#     chmod 755 install.sh &&\
#     bash install.sh &&\
#     rm install.sh &&\
#     yum install -y fuse.x86_64
# # Edit fuse config to allow root to access mount
# RUN echo "user_allow_other" >> /etc/fuse.conf


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
ENV WORKDIR=/home/serratus
CMD cd $WORKDIR

# Create serratus alias
RUN echo 'alias serratus="sudo -u serratus"' >> /root/.bashrc
RUN echo 'alias fusera="sudo -u serratus /usr/local/bin/fusera"' >> /root/.bashrc

## SRA Mounting Directory
#  First create a directory to mount data.
#  does not work
#  RUN serratus mkdir ./sra-mount

# # Testing ---

# # SRA toolkit
# vdb-config -i & read -t 3 ;  kill $!

# # Create aws key (read from ec2 IAM) (run at start-up of container)
# AWS_ACCESSKEYID=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.AccessKeyId' )
# AWS_SECRETKEY=$(curl http://169.254.169.254/latest/meta-data/identity-credentials/ec2/security-credentials/ec2-instance/ | jq -r '.SecretAccessKey' )
# echo User name,Password,Access key ID,Secret access key,Console login link > key.csv
# echo default,,$AWS_ACCESSKEYID,$AWS_SECRETKEY, >> $WORKDIR/key.csv

# vdb-config --accept-aws-charges yes \
#   --report-cloud-identity yes \
#   --set-aws-credentials $WORKDIR/key.csv

# fastq-dump SRR11183200


# # Example mount command:

# ## Command works; unsure if it's an S3-mount or goes via network
# ## need to confirm!
# wget https://serratus-public.s3-us-west-2.amazonaws.com/tmp/prj_phs710EA_test.ngc
# wget https://serratus-public.s3-us-west-2.amazonaws.com/tmp/SRR_Acc_List_Example.txt

# fusera mount -t ./prj_phs710EA_test.ngc -a ./SRR_Acc_List_Example.txt  ./sra-mount > output.log 2>&1 &

# fusera mount -t ./prj_phs710EA_test.ngc -a SRR5572888 ./sra-mount > output.log 2>&1 &

# serratus mkdir sra-mount

# fusera mount \
#  --a SRR5754947 \
#   ./sra-mount 2>&1 &

#  --aws-profile Default \
# --location s3.us-east-1 \
# -v \
# #Note: still requires switch to bash to even see the file
# # cannot read file yet, only see it exists (maybe)

# ## Unmount
# fusera unmount ./sra-mount


# sudo -u serratus /usr/local/bin/fusera mount -a "" ./sra-mount > output.log 2>&1 & 

# # sudo -u serratus /usr/local/bin/fusera mount -a "https://sra-pub-src-1.us-east-1.s3.amazonaws.com/SRR8388666" ./sra-mount > output.log 2>&1 & 





# #  Next run this fusera command to mount runs using the included repository key and accession list.

# fusera mount -t ./prj_phs710EA_test.ngc -a ./SRR_Acc_List_Example.txt ./example_mount > output.log 2>&1 &

# serratus mkdir $WORKDIR/sra

#   The example_mount directory will contain four SRR directories containing source data files and output.log will contain any messages or errors from the mount process.


# ## READING DATA

#   The mounted files can be read by commons tools as normal

# samtools view -H ~/example/example_mount/SRR1270860/NA19451.unmapped.ILLUMINA.bwa.LWK.low_coverage.20120522.bam


# ## UNMOUNTING DATA
#   Run the following command to unmount the data files

# fusera unmount ~/example/example_mount
