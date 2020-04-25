# Contributing to Serratus

  * [Open Science Collaboration](#open-science-collaboration)
    + [Contact the team](#contact-the-team)

  * [Bioinformatics pipeline](#bioinformatics-pipeline)
  * [AWS architecture](#aws-architecture)
  * [Repository organization](#repository-organization)

  * [Getting Started](#getting-started)
    + [Finding a Task](#finding-a-task)
    + [Creating a new Task](#creating-a-new-task)
    + [Using git](#using-git)
      - [Development is on Branches](#development-is-on-branches)
      - [Pull Request](#pull-request)

  * [Running an Experiment](#running-an-experiment)
    + [Experiment organization](#experiment-organization)
    + [Experiment Template](#experiment-template)
    + [Depositing Data](#depositing-data)

  * [Serratus Data](#serratus-data)
    + [`~/notebook` : Experiment associated data](#---notebook----experiment-associated-data)
    + [`~/out` : Serratus alignment output](#---out----serratus-alignment-output)
    + [`~/seq` : Coronavirus pan-genomes](#---seq----coronavirus-pan-genomes)
    + [`~/sra` : SraRunInfo Tables (.csv.gz)](#---sra----sraruninfo-tables--csvgz-)
    + [`~/test-data` : example data for development](#---test-data----example-data-for-development)
    + [`~/var/` : Assorted nuts and bolts](#---var-----assorted-nuts-and-bolts)

  * [Production Containers and Code](#production-containers-and-code)
    + [Building containers](#building-containers)
    + [Uploading container images to AWS ECR](#uploading-container-images-to-aws-ecr)
    + [Run interactive an serratus-dl](#run-interactive-an-serratus-dl)
    + [Testing scheduler](#testing-scheduler)

  * [Data Release Policy](#data-release-policy)


## Open Science Collaboration
`Serratus` is an Open-Science project. Our aim is to create a 100% reproducible study with 100% transparent and freely available data.

_We welcome all scientists and developers to contribute._

### Contact the team
Email (ababaian AT bccrc DOT ca) or join our [Slack (type `/join #serratus`)](https://join.slack.com/t/hackseq-rna/shared_invite/zt-dwdg5uw0-TTcfrFagariqKpOSU_d6wg)

## Bioinformatics pipeline
< Workflow is under active development. >

## AWS architecture
![serratus-overview](img/serratus_overview.png)

The workhorse for `serratus` is currently the `C5.large` EC2 instance on AWS. Each node has 2 vCPU and 4 GB of memory. Every instance is a blank slate running `amazon linux 2` with `docker`. Workflow is encapsulated in a `container`. This allows for rapid and cheap scaling of the cluster.

---

## Repository organization

```
./serratus/
├── bin                     # Binaries and executable tools/modules         
├── containers              # Container make files and scripts for production
├── data                    # For _local_ storage of data (.gitignore)
├── doc                     # Documentation files
├── img                     # Visual assets and workflow diagrams
├── local                   # For _local_ storage of files (.gitignore)
├── notebook                # Shared electronic lab-notebook and associated data files.
├── packer                  # Creating standardized node images (AMI)
├── src                     # Source code for modules/tools used in Serratus
├── terraform               # Cloud resource definitions for cluster
├── CONTRIBUTING.md         
├── LICENSE
└── README.md
```

---


## Getting Started

Serratus requires [`git`](https://guides.github.com/introduction/git-handbook/) for version control and to manage contributions from many authors.

To download or clone the `serratus` repository locally use:
```
git clone https://github.com/ababaian/serratus.git
```

### Finding a Task

Development of `serratus` is ongoing. To find and solve an open development problem see our [Project Page](https://github.com/ababaian/serratus/projects/1). This is a prioritized list of "Open Tasks" that need to be done, "Tasks in Progress" which are currently being worked on by others, "Code Review" and "Completed Tasks".

Also you can browse all tasks which are organized as ["Issues" on github](https://github.com/ababaian/serratus/issues?q=).

Feel free to comment on any issue, even those you're not assigned to if you have a helpful suggestion.

If you'd like to work on a given task, simply add a comment saying this to the issue and it will be "Assigned" to you.

### Creating a new Task

If you have an idea you'd like to develop, would like to run an experiment or require additional documentation, let the other developers know what you're doing by creating an "Issue" on github. The general template to include initially is:

```
### Problem / Objective
< Briefly outline a problem you are solving / the research objectives and hypothesis you are testing >

### Proposed Solution / Methods
< How are you planning on solving the problem / experimental design to test the hypothesis >

### Additional Resources
< Outline any additional information you require to do this task or resources you'll need access to >

```

### Using git

#### Development is on Branches
In `git` there are independent working copies of the code-base called `branches`. Each contributor works on his/her own branch, and this does not interact with other branches.

The main or `master-branch` is for production use. This is the operational code and completed experiments. You should not add your changes to `master` as this can break serratus.

You should work on your own development branch. By convention name the branch `<feature_your_developing>-dev` or `<name>-dev`.

```
# Create a development branch called `sally-dev`
git branch sally-dev
# Move to your branch
git checkout sally-dev
```

Now you can modify any of the files in your local `~/serratus/`, when you're ready to save your work (say I modified `serratus/docker/Dockerfile`)

```
# Stage your changes
git add docker/Dockerfile

# Check the status of what changes are staged and what is not staged
git status

# Commit your changes to version control
git commit -m "Short description of what these changes accomplished"

# End of the work-day, back-up your changes to the online github repository
git pull # download others changes from the online repo
git push # upload your changes to the online repo 
```

#### Pull Request

When you've finished a larger unit of work, like completing a feature. Your `branch` changes can be combined with the `master` branch and into production code. This is called issuing a "Pull Request".

---

## Running an Experiment

It's important that every experimental result in `serratus` be 100% reproducible and all the data is freely available amongst contributors.

All results are documented in a shared electronic notebook: `serratus/notebook`

To assist in this, we suggest using `Jupyter Notebook` to document all the code in running an experiment. Alternatively you can use `markdown` or other ways where code and output can be distinguished.

### Experiment organization

Notebook entries naming and file convention. (`YYMMDD` is start date)
```
# Notebook Entry
serratus/notebook/YYMMDD_experiment_name.ipyb

# Scripts, small data files (>250 kb) and plots from experiment
serratus/notebook/YYMMDD/<filename>.Rscript
serratus/notebook/YYMMDD/<filename>.png
serratus/notebook/YYMMDD/<filename>.csv
...

# Large data files are stored on S3
s3://serratus-public/notebook/YYMMDD/bam/aligned.bam
s3://serratus-public/notebook/YYMMDD/large_data.Rdata
...
```

### Experiment Template
Copy the experiment template to start your experiment.
```
cp notebook/200401_template.ipynb notebook/200420_my_experiment_title.ipynb
```

### Depositing Data

Please contact us to be granted permission to deposit data 

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

### Building containers

```
# Start docker service on amazon linux 2
sudo yum install -y docker
sudo yum install -y git
sudo service docker start
```

```
# Download latest serratus repo
git clone https://github.com/ababaian/serratus.git; cd serratus/containers

# If you want to upload containers to your repository, include this.
export DOCKERHUB_USER='serratusbio' # optional
sudo docker login # optional

# Build all containers and upload them docker hub repo (if available)
./build_containers.sh

```

### Uploading container images to AWS ECR
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

### Run interactive an serratus-dl
```
sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:latest
```

### Testing scheduler
```
cd serratus/scheduler
sudo docker build -t scheduler:0 .
docker run -d --rm -p8000:8000 --name sch scheduler:0
curl -T /path/to/SraRunInfo.csv localhost:8000/add_sra_run_info
```

---

## Data Release Policy
To achieve our objective of providing high quality CoV sequence data to the global research effort, Serratus ensures:
- All software development is open-source and freely available (GPLv3)
- All sequencing data generated, raw and processed, will be freely available in the public domain in accordance with the [Bermuda Principles](https://en.wikipedia.org/wiki/Bermuda_Principles) of the Human Genome Project.