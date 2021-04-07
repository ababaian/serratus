# serratus_genbankParser.Rmd

## Usage

not formatted as script yet
    
### Dependecies
    Modules: R 4.0, BiocManager(devel), genbankr(devel), taxize, data.table, rlist, devtools, GenomicRanges, Biostrings, tidyverse
    Files: cov0.gb, cov0.duplicates and cov0.id99.uc.

## Description

The genbankParser is designed to be run as a standalone script to generate a formatted and cleaned csv table of the covid pan-genome from genbank input. It deals with most hostTaxonId mapping errors, and attempts to infer these hostTaxonIds for duplicate and highly homologous entries by checking if clusters/duplicate all provide the same hostTaxonId and if so inferring it for those where none was prvided (in a new column). 

Three meta-data files are required: the covid pan-genome from genbank [file](https://serratus-public.s3.amazonaws.com/seq/cov0/cov0.gb), the list of duplicates (generated here [file](https://github.com/ababaian/serratus/blob/master/notebook/200420_cov2_pangenome.ipynb) and a table containing homology information [file](https://serratus-public.s3.amazonaws.com/seq/cov2r/cov0.id99.uc) Filepaths have to be edited into the code as it stands
