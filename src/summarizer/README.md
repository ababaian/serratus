# serratus_summarizer.py

## Usage

python3 serratus_summarizer.py InputFileName SummaryFileName OutputSAMFileName [OutputTripletFileName]
    
### Positional arguments:
    1. InputFileName      Input filename, SAM filename or /dev/stdin
    2. SummaryFileName    Summary filename, text format
    3. OutputSAMFileName     Output filename, usually /dev/stdout or - for none. Input is echo'd to this file.
    
### Optional positional argument:
    4. TipletFileName     Output filename to store tsv file with (Label, Start, Length) for each alignment.
    
### Dependecies
    Modules: None (stand-alone python3).
    Files: acc_len_taxid.txt and taxid_desc.txt.

### Usage in pipeline

    $ bowtie2 <args> \
        | python3 serratus_summarizer.py /dev/stdin summary.txt /dev/stdout tiplets.txt \
        | samtools <args>

### Stand-alone usage

    $ python3 serratus_summarizer.py SRR1234.sam summary.txt - # notice trailing "-"

## Description

The summarizer is designed to be inserted into the mapping pipeline to generate a summary of hits when a dataset completes. It acts like the Linux 'tee' command: input is forwarded to output while collecting data.

Two meta-data files are required: acc_len_taxid.txt and taxid_desc.txt which are derived from NCBI data (Genbank records for pan-genome accessions and the Taxonomy database, respectively). By default, they are expected in the current directory. This can be overridden by exporting the SUMZER_DIR environment variable. These files can be obtained from here:

    s3://serratus-public/var/acc_len_taxid.txt
    s3://serratus-public/var/taxid_desc.txt

