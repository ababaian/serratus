# serratus
![Serratus Mountain in Squamish, BC. Canada](img/serratus_logo.png)

*Serratus Mountain, Squamish,BC*

### Background
The SARS-CoV-2 pandemic will infect millions and has already crippled the global economy. 

While there is an intense research effort to sequence SARS-CoV-2 isolates to understand the evolution of the virus in real-time, our understanding of where it originated is limited by the sparse characterization of other members of the Coronaviridae family (only 53/436 CoV sp. Genomes are available).
 
We are re-analyzing all RNA-sequencing data in the NCBI Short Read Archive to discover new members of Coronaviridae. Our initial focus is mammalian RNA-sequencing libraries followed by avian/vertebrate, metagenomic, and finally all 1.12M entries (5.72 petabytes).

### Architecture
![serratus-overview](img/serratus_overview.png)


## Repository organization

`data/README.md`: README.md outlining location/acquisition of serratus data

`docker/`: Container make files and job scripts

`img/`: Architecture/workflow diagrams

`packer/`: Standardized node images (ami)

`scheduler/`: Code for `serratus` head-node and sraRunInfo management

`terraform/`: Cloud resources / pipeline management


## Useful links
- **S3 Bucket:** s3://serratus-public/ (public-readable)

#### AWS-specific
- [AWS Batch workflow - Introduction](https://aws.amazon.com/blogs/compute/building-high-throughput-genomics-batch-workflows-on-aws-introduction-part-1-of-4/)
- [AWS Batch workflow - github page](https://github.com/aws-samples/aws-batch-genomics)
- [SRAtoolkit in Cloud Computing](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/)
- [NCBI SRA Data on S3](https://registry.opendata.aws/ncbi-sra/)
- [S3 transfer optimization](https://docs.aws.amazon.com/cli/latest/topic/s3-config.html)
- [Paper on analyzing EC2 costs (2011)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0026624)
- [Pushing the limits of Amazon S3 Upload Performance](https://improve.dk/pushing-the-limits-of-amazon-s3-upload-performance/)
- [Clever SRA alignment pipeline](https://github.com/FredHutch/sra-pipeline
)

#### SARS-CoV-2
- [SARS-CoV-2 UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=wuhCor1)
- [Interpretable detection of novel human viruses from genome sequencing data](https://www.biorxiv.org/content/10.1101/2020.01.29.925354v3.full.pdf)
- [Virus detection from RNA-seq: proof of concept](https://www.ncbi.nlm.nih.gov/pubmed/21603639)

#### Bloom Filters
- [Bigsi: Bloom filter indexing of SRA/ENA for organism search](https://github.com/phelimb/bigsi)
- [Fast Search of Thousands of Short-Read Sequencing Experiments](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4804353/)
- [Ultra-fast search of all deposited bacterial and viral genomic data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6420049/)


# Setting up and running Serratus

### 0) Dependencies

#### AWS account

1. Sign up for an AWS account (you can use the free tier)
2. [Create an IAM Admin User with Access Key](https://docs.aws.amazon.com/IAM/latest/UserGuide/getting-started_create-admin-group.html). For **Access type**, use **Progammatic access**.
3. Note the Access Key ID and Secret values.
4. Create a [EC2 keypair](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair) in `us-east-1` region. Retain the name of the keypair and the `.pem` file. Configure your `ssh` for easy AWS access(change `serratus.pem` to your identity file).

`~/.ssh/config`: Add these lines
```
Host *.compute.amazonaws.com *.compute-1.amazonaws.com aws_*
     User ec2-user
     IdentityFile ~/.ssh/serratus.pem
     StrictHostKeyChecking no
     UserKnownHostsFile /dev/null
```

#### Packer

1. [Download Packer](https://packer.io/downloads.html) as a binary. Extract it to a PATH directory (`~/.local/bin`)

#### Terraform

1. [Download Teraform](https://www.terraform.io/downloads.html) (>= v0.12.24) as a binary. Extract it to a PATH directory (`~/.local/bin`)

### 1) Build Serratus AMIs with Packer

Pass AWS credentials to pipeline via environmental variables
```
export AWS_ACCESS_KEY_ID="your_access_key"
export AWS_SECRET_ACCESS_KEY="your_secret_key"
```

Use packer to build the serratus instance image (AMI)
```
cd serratus/packer
/path/to/packer build docker-ami.json
cd ../..
```

This will start up a t3.nano, build the AMI, and then terminate it.  Currently this takes about 2 minutes, which should cost well under a penny. The final line of STDOUT will be the region and AMI. Retain this information

Current stable AMI: `us-east-1: ami-04c1625cf0bcb4159`

###  2) Build Serratus resources with Terraform

#### Set Terraform variables

Open `terraform/main/terraform.tfvars` in a text editor. Set these variables
 * `dev_cidrs`: Your public IP, followed by "/32". Use: `curl ipecho.net/plain; echo`
 * `key_name`: Your EC2 key pair name
 * `dockerhub_account`: (optional). Change this to your docker hub account to build your own images. Default images are in `serratusbio` organization.


#### Create Serratus resources

Navigate to the top-level module and run `terraform` initialization and apply. Retain the scheduler DNS address (last output line).

```
cd terraform/main
terraform init
terrafform apply
cd ../..
```
At the time of writing, this will create:

  * a t3.nano, for the scheduler, with an Elastic IP
  * an S3 bucket, to store intermediates
  * an ASG for serratus-dl, using c5.large with 50GB of gp2.
  * An ASG for serratus-align, using c5.large
  * An ASG for serratus-merge, using t3.small
  * Security groups and IAM roles to tie it all together.

All ASGs have a max size of 1.  This can all be reconfigured in terraform/main/main.tf.

At the end of `tf apply`, it will output the scheduler's DNS address.  Keep this for later.

### 3) Open SSH tunnel to the scheduler

The scheduler exposes ports 3000/8000/9090.  This port is *not* exposed to the public internet. You will need to create an SSH tunnel to allow your local web-browser and terminal to connect.  

```
./create_tunnel.sh
```
Open a web browser for UI: [Status Page: http://localhost:8000/jobs/](http://localhost:8000/jobs/) [Grafana: http://localhost:3000/jobs/](http://localhost:8000/jobs/) [http://localhost:8000/jobs/](http://localhost:3000) [Prometheus: http://localhost:8000/jobs/](http://localhost:9090)

May take a few minutes to boot. Make tea.

### 5) Loading SRA Accesions into Serratus

Once the scheduler is online, you can curl SRA accession data in the form of a `SraRunInfo.csv` file (NCBI SRA > `Send to: File`).

```
curl -s -X POST -T /path/to/SraRunInfo.csv localhost:8000/jobs/add_sra_run_info/
```

This should respond with a short JSON indicating the number of rows inserted, and the total number in the scheduler.

In your web browser, refresh the status page.  You should now see a list of accessions by state. If ASGs are online, they should start processing immediately.  In a few seconds, the first entry will switch to "splitting" state, which means it's working.

### 6) Launch cluster nodes (currently manual)

With data loaded into the scheduler, manually set the number of `serratus-dl` & `serratus-align` nodes to process the data.



# Data Release Policy
To achieve our objective of providing high quality CoV sequence data to the global research effort, Serratus ensures:
- All software development is open-source and freely available (GPLv3)
- All sequencing data generated, raw and processed, will be freely available in the public domain in accordance with the [Bermuda Principles](https://en.wikipedia.org/wiki/Bermuda_Principles) of the Human Genome Project.

# Contributing
`Serratus` is an Open-Science project. We welcome all scientists to contribute.
Email (ababaian AT bccrc DOT ca) or join [Slack](https://join.slack.com/t/hackseq-rna/shared_invite/zt-cs2f4dy8-QZ92T8E1O_Lwrse18yXWEA)