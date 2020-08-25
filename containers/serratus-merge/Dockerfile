FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='merge'
ARG VERSION='0.1.4'

# Software Information
# ENV SAMTOOLSVERSION='1.10'  # from serratus-base


#==========================================================
# Dependencies ============================================
#==========================================================

#==========================================================
# Install Software ========================================
#==========================================================

#RUN curl https://serratus-public.s3.amazonaws.com/var/acc_len_taxid.txt > \
#  ./acc_len_taxid.txt
#RUN curl https://serratus-public.s3.amazonaws.com/var/taxid_desc.txt > \
#  ./taxid_desc.txt

FROM amazonlinux:2 AS runtime

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: merge container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools"

# Environment setup
ENV BASEDIR=/home/serratus

# aws cli, plus dependencies
# -merge has its own python dependency, so do a full python/pip install
RUN yum -y install python3 perl ps jq util-linux \
 && alias python=python3 \
 && curl -O https://bootstrap.pypa.io/get-pip.py \
 && python3 get-pip.py \
 && rm get-pip.py \
 && pip install boto3 awscli \
 && yum clean all \
#   aws configuration
 && aws configure set default.s3.multipart_threshold 4GB \
 && aws configure set default.s3.multipart_chunksize 4GB

# samtools
COPY --from=build_base /usr/local/bin/samtools /usr/local/bin/

#==========================================================
# Serratus Initialize =====================================
#==========================================================
WORKDIR /home/serratus
COPY worker.sh ./
COPY serratus-merge/*sh ./
COPY serratus-merge/*py ./

# sumbler exe
#RUN aws s3 cp s3://serratus-public/sumbler/sumbler ./ \
# && chmod 755 sumbler
RUN curl https://serratus-public.s3.amazonaws.com/var/sumbler > \
  ./sumbler \
 && chmod 755 sumbler

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "merge", "/home/serratus/serratus-merge.sh"]
