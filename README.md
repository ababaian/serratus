# serratus
![serratus](img/serratus_logo.png)

AWS framework for cheap and fault-tolerant DNA/RNAs-seq alignment.

## Repo folders

**/ami**: Make scripts for instance images (ami) and their respective ami-id

**/img**: Diagrams of serratus workflows

**/scripts**: Defined units of work (jobs) performed in serratus

# Project S3
**s3://serratus-public/**

**$S3/resources**: Genome indices
- hgr1-bt2/ : human rDNA bowtie2 indexed

**$S3/example-data**: Toy data for testing
- bam/ : aligned bam files for breaking into blocks
- fq/  : sequencing reads of various length


## Useful links
- **S3 Bucket:** s3://serratus-public/ (public-readable)
- [Alpine Linux AMI](https://github.com/mcrute/alpine-ec2-ami)
- [AWS Batch workflow - Introduction](https://aws.amazon.com/blogs/compute/building-high-throughput-genomics-batch-workflows-on-aws-introduction-part-1-of-4/)
- [AWS Batch workflow - github](https://github.com/aws-samples/aws-batch-genomics)
- [FSx Shared file-system for HPC](https://aws.amazon.com/blogs/storage/using-amazon-fsx-for-lustre-for-genomics-workflows-on-aws/)
- [SRA in the Cloud -- Use SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/)
- [SRA Data Registry on S3](https://registry.opendata.aws/ncbi-sra/)
- [S3 transfer optimization](https://docs.aws.amazon.com/cli/latest/topic/s3-config.html)
- [Paper on analyzing EC2 costs (2011)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0026624)
- [Pushing the limits of Amazon S3 Upload Performance](https://improve.dk/pushing-the-limits-of-amazon-s3-upload-performance/)
- [Clever SRA alignment pipeline](https://github.com/FredHutch/sra-pipeline
)
- [Interpretable detection of novel human viruses from genome sequencing data](https://www.biorxiv.org/content/10.1101/2020.01.29.925354v3.full.pdf)

#### Architecture

![serratus-overview](img/serratus_overview.png)

# Getting up and running

## Building AMIs with Packer

First, [download Packer](https://packer.io/downloads.html).  It comes as a single
binary which you can just unzip.  I extracted it to `~/.local/bin` so that it ended
up on my PATH.

Next, use it to build the AMI: `/path/to/packer build serratus/packer/docker-ami.json`

This will start up a t3.nano, build the AMI, and then terminate it.  Currently this
takes about 2 minutes, which should cost well under a penny.

## Getting started with Terraform

The top-level module is in serratus/terraform/cluster.  If you're OK with the
defaults, you shouldn't need to change anything.  Just cd there, run `tf init`,
and then `tf apply` to create some infrastructure.

