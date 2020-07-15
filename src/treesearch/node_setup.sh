#!/bin/bash

set -e

. "/etc/parallelcluster/cfnconfig"

MASTER=1

# differentiate master from compute node
case "${cfn_node_type}" in
  MasterServer)
    MASTER=1
  ;;
  ComputeFleet)
    MASTER=0
  ;;
  *)
  ;;
esac

# install basics
# yum update
yum -y install wget cmake cmake3 zlib-devel gzip unzip flex bison

# slave cmake3 as the defaul
alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake 10 \
--slave /usr/local/bin/ctest ctest /usr/bin/ctest \
--slave /usr/local/bin/cpack cpack /usr/bin/cpack \
--slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake \
--family cmake
alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake3 20 \
--slave /usr/local/bin/ctest ctest /usr/bin/ctest3 \
--slave /usr/local/bin/cpack cpack /usr/bin/cpack3 \
--slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake3 \
--family cmake

# install usearch
wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz -O usearch.gz
gunzip usearch.gz
chmod ugo+x usearch
mv usearch /usr/local/bin

if [[ $MASTER -eq 1 ]]
then
  # safety copy of slurm conf
  cp /opt/slurm/etc/slurm.conf /opt/slurm/etc/slurm.conf.backup
  # enable accounting via file in slurm
  # WARNING probably change all this if you decide to use a DB for accounting instead
  sed -i -E 's/^JobCompType=.*$/JobCompType=jobcomp\/filetxt/g' /opt/slurm/etc/slurm.conf
  sed -i -E 's/^#JobCompLoc=/JobCompLoc=\/shared\/slurm\/job_completions/g' /opt/slurm/etc/slurm.conf
  sed -i -E 's/^#JobAcctGatherType=.*$/JobAcctGatherType=jobacct_gather\/linux/g' /opt/slurm/etc/slurm.conf
  sed -i -E 's/^#AccountingStorageType=.*$/AccountingStorageType=accounting_storage\/filetxt/g' /opt/slurm/etc/slurm.conf
  sed -i -E 's/^#AccountingStorageLoc=.*$/AccountingStorageLoc=\/shared\/slurm\/accounting/g' /opt/slurm/etc/slurm.conf

  # set up the accounting files
  mkdir -p /shared/slurm
  touch /shared/slurm/accounting
  chown slurm:slurm /shared/slurm/accounting
  touch /shared/slurm/job_completions
  chown slurm:slurm /shared/slurm/job_completions

  cd /home/centos
  # install miniconda
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  chown centos:centos miniconda.sh
  chmod u+x miniconda.sh
  runuser -l  centos -c './miniconda.sh -b'
  runuser -l  centos -c './miniconda3/condabin/conda init'
  rm miniconda.sh

  # install pargenes
  git clone --recursive https://github.com/BenoitMorel/ParGenes.git
  chown -R centos:centos ParGenes/
  cd ParGenes
  git checkout tags/v1.1.2
  ./install.sh 18
  cd -

  # install snakemake
  runuser -l  centos -c './miniconda3/condabin/conda install -y -c conda-forge -c bioconda snakemake'

  # get the pipeline
  git clone --recursive https://github.com/lczech/nidhoggr.git
  chown -R centos:centos nidhoggr/

fi