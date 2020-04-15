FROM serratus-base:latest AS build_base

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Container Information
ARG TYPE='merge'
ARG VERSION='0.1.2'

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
COPY docker/worker.sh ./
COPY docker/serratus-merge/*sh ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["./worker.sh", "merge", "/home/serratus/serratus-merge.sh"]
