# serratus

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


#### Architecture

![serratus-overview](img/serratus_overview.png)