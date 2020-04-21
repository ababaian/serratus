# Contributing to Serratus

## Bioinformatics pipeline
< Placeholder for bioinformatics description >

## AWS architecture
< Placeholder for architecture description >

## Open Science Collaboration
< Placeholder to outline how data/experiments are shared >

---

## Repository organization

- `data/`: For _local_ storage of data.
- `docker/`: Container make files and scripts for production `serratus`.
(will be renamed to `containers`)
- `img/`: Visual assets and workflow diagrams
- `packer/`: Creating standardized node images (ami)
- `scheduler/`: Code for `serratus` head-node and sraRunInfo management.
(will be moved to `containers`)
- `terraform/`: Cloud resource definitions for creating/running cluster.
- `notebook/`: Shared electronic lab-notebook entries and associated data files.
- `runs/`: Scripts for analyzing the summary outputs of `serratus`

---

## Serratus Data

**Serratus Bucket(`~`)**: `s3://serratus-public/`

All data files are stored on our AWS S3 bucket. To access:

1) Using `aws-cli` : `aws s3 cp s3://serratus-public/<file_path> ./`

2) Using `wget` : `https://serratus-public.s3-us-east-1.amazonaws.com/<file_path> ./`

To browse directory: `aws s3 ls s3://serratus-public/`

### `~/notebook` : Experiment associated data
For each electronic lab notebook entry, data associated with that run can be stored in this directory. Each folder is a date (`YYMMDD`) corresponding to the date of the notebook file. For example

The data for the experiment (`serratus/notebook/200411_CoV_Divergence_Simulations.ipynb`](https://github.com/ababaian/serratus/blob/master/notebook/) is found in `s3://serratus-public/notebook/200411/`.


### `~/out` : Serratus alignment output
Currently output from ~50 librarires using standard alignment
- `~/out/bam` : Aligned output file + index, SRA accession named
- `~/out/flagstat` : Flagstat files for bam files

### `~/seq` : Coronavirus pan-genomes
Reference genome sets and their associated index files.

- `~/seq/cov0`  : All CoV sequences from NCBI
	- NCBI search: `"(Coronaviridae) AND "viruses"[porgn:txid10239]"`
	- Date Accessed: 2020/03/30
	- Results: 33296

- `~/seq/cov1r`  : Initial pan-coronavirus genome
	- Based off of `cov0` with non-CoV accessions removed and polyA[10+] masked
	- See: `~/serratus/notebook/200408_cov1_pangenome.ipynb` for make commands
	- Date: 2020/04/08
	- `cov01r` contains reverse non-compliment control sequences

- `~/seq/cov2r`  : Refined pan-coronavirus genome
	- Based off of `cov0`
	- Removed poly-nt tracts of 10+
	- Blacklisted 6 non-CoV accessions
	- Pruned 
	- See: `~/serratus/notebook/200420_cov2_pangenome.ipynb` notebook for commands
	- Date: 2020/04/20
	- `cov2.fa`  : Masked pan-genome 
	- `cov2r.fa` : Masked pan-genome with reverse non-compliment controls

- `~/seq/hgr1`  : Human rDNA testing sequence
	- From [this publication](https://www.biorxiv.org/content/10.1101/118760v2)


### `~/sra` : SraRunInfo Tables (.csv.gz)
SRA Accession and Run Information master tables. Accessed via SRA website and the following basic filter:

```
"type_rnaseq"[Filter] AND cluster_public[prop] AND "platform illumina"[Properties] AND "cloud s3"[Properties] NOT "scRNA"[All Fields] AND <SUBFILTER>
```

- **Test Data Set**
	- Mammals and CoV+ swabs for testing pipeline
	- SARS-CoV-2: `PRJNA616446`
	- Felis catus: `PRJNA432069`
	- Homo sapiens (HCT116): `PRJEB29794`
	- Macaca fascicularis: `PRJNA553361`
	- Mus musculus: `PRJNA553361`
	- Date Accessed: 2020/04/07
	- Results: 49 libraries

- **Non-Human, Non-Mouse Mammals**
	- `BASE AND "Mammalia"[Organism] NOT "Homo sapiens"[Organism]) NOT "Mus musculus"[orgn]`
	- Date Accessed: 2020/03/28
	- Results: 66926, 0.15 PB

- **Human**
	- `BASE AND "Homo sapiens"[Organism]`
	- Date Accessed: 2020/03/05
	- Results: 520257, 4.75 PB

- **Mouse**
	- `BASE AND "Mus musculus"[orgn]`
	- Results: 539233
	- Not accessed

- **Vertebrates, Non-mammal**
	- `BASE NOT "Mammalia"[Organism] NOT "Homo sapiens"[Organism] NOT "Mus musculus"[orgn]`
	- Date Accessed: 2020/03/29
	- Results: 74532, 0.115 PB

- **Invertebrates**
	- `BASE NOT "Vertebrata"[Organism]`
	- Date Accessed: 2020/03/30
	- Results: 403639, 0.7 PB

- **HCT116 RNAseq**
	- For testing; ca. 1000 entries of human HCT116 cell line

### `~/test-data` : example data for development

**Sequence Files**
- `../bam/` : aligned bam files for breaking into blocks
- `../bam-block` : bam file output of fq-blocks requiring merging
- `../fq/`  : sequencing reads of various length
- `../fq-block` : fq files broken into 'blocks'
- `../out` : Example output data of re-aligned reads

### `~/var/` : Assorted nuts and bolts

---

## Production Containers and Code

```
## Container Image Hierarchy
# serratus-base
# |    |--serratus-dl
# |    |--serratus-align
# |    |     |--serratus-align-hgr1
# |    |     |--...
# |    |--serratus-merge
# serratus-scheduler
# serratus-grafana
# serratus-prometheus
```

#### Building containers

```
# Start docker service on amazon linux 2
sudo yum install -y docker
sudo yum install -y git
sudo service docker start
```

```
# Download latest serratus repo
git clone https://github.com/ababaian/serratus.git; cd serratus

# If you want to upload containers to your repository, include this.
export DOCKERHUB_USER='serratusbio' # optional
sudo docker login # optional

# Build all containers and upload them docker hub repo (if available)
./build.sh

```

#### Uploading container images to AWS ECR
Paste resulting command in terminal to authenticate (not implemented)

```
## NOTE: Depricated. Jeff builds to dockerhub directly
##       these instructions are for 
#aws ecr get-login --region us-east-1 --no-include-email

## Tag images with ECR information
# sudo docker tag serratus-base:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker tag serratus-dl:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker tag serratus-align:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker tag serratus-merge:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge

## Push images to ECR (rather host on dockerhub)
## Make sure to initialize each repository on your AWS ECR account
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge
```

#### Run interactive an serratus-dl
```
sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:latest
```

#### Testing scheduler
```
cd serratus/scheduler
sudo docker build -t scheduler:0 .
docker run -d --rm -p8000:8000 --name sch scheduler:0
curl -T /path/to/SraRunInfo.csv localhost:8000/add_sra_run_info
```

---

## Contributing
### Using Git/Github
< Best practices for using `git` >

### Finding an open Task

### Creating a new Task

### Running an experiment

### Accessing Data

### Depositing Data

### Report a bug / Suggest Feature

