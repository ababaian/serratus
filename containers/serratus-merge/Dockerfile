FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='merge'
ARG VERSION='0.1.4'

# Software Information
# ENV SAMTOOLSVERSION='1.10'  # from serratus-base

# Additional Metadata
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="serratus: alignment container"
LABEL software.license="GPLv3"
LABEL tags="aws-cli, samtools"

#==========================================================
# Dependencies ============================================
#==========================================================

#==========================================================
# Install Software ========================================
#==========================================================

#==========================================================
# Serratus Initialize =====================================
#==========================================================
WORKDIR /home/serratus
COPY worker.sh ./
COPY serratus-merge/*sh ./
COPY serratus-merge/*py ./
#RUN curl https://serratus-public.s3.amazonaws.com/var/acc_len_taxid.txt > \
#  ./acc_len_taxid.txt
#RUN curl https://serratus-public.s3.amazonaws.com/var/taxid_desc.txt > \
#  ./taxid_desc.txt

FROM amazonlinux:2 AS runtime

COPY --from=build_base / /

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "merge", "/home/serratus/serratus-merge.sh"]
