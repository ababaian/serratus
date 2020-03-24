# Docker Base: amazon linux2
FROM amazonlinux:2 AS build_base

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
COPY scripts/worker.sh ./
COPY scripts/serratus-merge/* ./

#==========================================================
# ENTRYPOINT ==============================================
#==========================================================
ENTRYPOINT ["/bin/bash", "/home/serratus/serratus-merge.sh"]