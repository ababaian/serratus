# Contributing to Serratus

**Serratus** is an Open-Science project. Our aim is to unlock all of Earth's viruses in public data, freely and with 100% transparency.

_We welcome all scientists and developers to contribute._

## Collaborating with Serratus
Our data is free and available for everyone. If you require assistance with accessing or using Serratus data, please let use know we would love to help. We are set up and capable to run assemblies and retrieve full-length viruses where available. We ask for nothing in return.

# Ongoing Project Roles

Below are a few outlines of things which need development. Skill-sets listed are all optional, the most important trait is to be a "self-starter" and passionate to learn. We're all volunteers and will help one another to get there.

### Command Line Interface Design

The `Serratus` command line interface is currently a mixture of scripts, jupityer notebooks and tweaking terraform configuration files. The goal is to create a stream-lined command-line interface for running `Serratus` and writing the associated documentation. The outcome will be a readily deployable workflow.

Skills: BASH scripting, written communication, AWS systems, Terraform deployment.

### Web Interface Design

Our web interface: www.serratus.io , is under active development. We are looking for new and intereting ways to display and make sense of this ocean of data. In particular we're developing an "RdRP" focused characterization tool called [`palmID`](https://serratus.io/palmid).

Skills: JavaScript, R, Docker, SQL, UX Design.

### Virus Host Interaction Networks

The `Serratus` data spans all of Earth, all domains of life, and thousands of ecological niches. We are trying to make sense of the global distribution of RNA viruses using a "Data-Driven" philosophy (unsupervised). The goal is to compile available virus-host data and use it for the computational annotation of novel sequences. You can get an idea of how this will work roughly in the [`palmID` web interface](https://serratus.io/palmid?hash=3xample)

Skills: SQL/Graph databases, Virus-host modeling, Ecoinformatics, R, Machine Learning.

### Virus-X Characterization

There are thounsands of novel viruses and virus families uncovered through our RdRP search. If you would like to focus on one particular group of viruses, characterize them and make sense of the data, we can help :)

Skills: Virology domain-knowledge, phylogenetics, bioinformatics, writing.

### Bioinfo-magician

At the heart of `Serratus` is finding ways to optimize sequence-search at a massive scale. If you're keen to get into the "guts" of a bunch of software and look for ways to speed things up, this is always a goal.

Skills: C/C++, sequence search, AWS systems, optimization, wizardry.

### ????

Everything we didn't think of. If you have a cool idea and think it can fit in, we're happy to hear it and work together.

# How to Contribute

Skim this document and [join our slack](https://join.slack.com/t/hackseq-rna/shared_invite/zt-ewlzh9qf-SiNkxvvTJflcutFN0h5jIQ) and say hello!

## Project Development Overview

- Development is [done through `git`](https://github.com/ababaian/serratus/wiki/Using-the-Repository).
- Data is hosted on [a public AWS S3 bucket](https://github.com/ababaian/serratus/wiki/Access-Data-Release).
- [Deployment of Serratus resources](https://github.com/ababaian/serratus/wiki/Running-Serratus-on-AWS) can be done on a personal AWS account.
- Experiments are [run using Jupyter Notebooks](https://github.com/ababaian/serratus/wiki/Running-an-Experiment).
- Contact the team: join our [Slack: (type `/join #serratus`)](https://join.slack.com/t/hackseq-rna/shared_invite/zt-dwdg5uw0-TTcfrFagariqKpOSU_d6wg) or email `artem AT rRNA DOT ca`

## Find an existing task

To find and solve an open development problem see our [Project Page](https://github.com/ababaian/serratus/projects/1). This is a prioritized list of "Open Tasks" that need to be done, "Tasks in Progress" which are currently being worked on by others, "Code Review" and "Completed Tasks".

Also you can browse all tasks which are organized as ["Issues" on github](https://github.com/ababaian/serratus/issues?q=).

Feel free to comment on any issue, even those you're not assigned to if you have a helpful suggestion.

If you'd like to work on a given task, simply add a comment saying this to the issue and it will be "Assigned" to you.

## Create a new task

If you have an idea you'd like to develop, would like to run an experiment or require additional documentation, let the other developers know what you're doing by [creating an "Issue" on Github](https://github.com/ababaian/serratus/issues/new). The general template to include initially is:

```
### Problem / Objective
< Briefly outline a problem you are solving / the research objectives and hypothesis you are testing >

### Proposed Solution / Methods
< How are you planning on solving the problem / experimental design to test the hypothesis >

### Additional Resources
< Outline any additional information you require to do this task or resources you'll need access to >

```

# Project policies

## Authorship guidelines

There is no formal structure to the Serratus team, everyone is encouraged to take full ownership of the components they develop. For publications authorship is determined by the [ICMJE Guidelines](http://www.icmje.org/recommendations/browse/roles-and-responsibilities/defining-the-role-of-authors-and-contributors.html). Specific emphasis is placed on that authors must be directly involved in the collaboration.

## Data and Software Policies

To achieve our objective of providing high quality CoV sequence data to the global research effort, Serratus ensures:

- All software development is open-source and freely available (GPLv3)
- We adopt the [INSDC Release Policy (2002)](http://www.insdc.org/policy.html) for Serratus data, databases and derivative analyses
- All data generated, raw and processed, will be freely available in the public domain in accordance with the [Bermuda Principles](https://en.wikipedia.org/wiki/Bermuda_Principles) of the Human Genome Project.
