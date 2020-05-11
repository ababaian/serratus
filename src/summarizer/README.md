# serratus_summarizer.py

## Usage

python3 serratus_summarizer.py InputFileName SummaryFileName OutputSAMFileName [OutputTripletFileName]
    
### Positional arguments:
    1. InputFileName      Input filename, SAM filename or /dev/stdin
    2. SummaryFileName    Summary filename, text format
    3. OutputSAMFileName     Output filename, usually /dev/stdout or - for none. Input is echo'd to this file.
    
### Optional positional argument:
    4. TinyhitFileName     Output filename to store tinyhits.
    
### Dependecies
    Modules: None (stand-alone python3).
    Files: acc_len_taxid.txt and taxid_desc.txt.

### Usage in pipeline

    $ bowtie2 <args> \
        | python3 serratus_summarizer.py /dev/stdin summary.txt /dev/stdout tiplets.txt \
        | samtools <args>

### Stand-alone usage

    $ python3 serratus_summarizer.py SRR1234.sam summary.txt - th.tsv # notice "-" for the output filename

## Description

The summarizer is designed to be inserted into the mapping pipeline to generate a summary of hits when a dataset completes. It acts like the Linux 'tee' command: input is forwarded to output while collecting data.

Two meta-data files are required: acc_len_taxid.txt and taxid_desc.txt which are derived from NCBI data (Genbank records for pan-genome accessions and the Taxonomy database, respectively). By default, they are expected in the current directory. This can be overridden by exporting the SUMZER_DIR environment variable.

Tinyhit is a new output format designed to be small as possible while enabling summary reports to be generated or re-generated. We expect the summarizer to evolve rapidly, and it would be better to work from a more compact format if possible. There are five fields: 1=Label, 2=Start, 3=Length, 4=diffs, 5=ee. Diffs is edit distance (NM:i tag in SAM record), ee=number of expected errors in the alignment calculated from the Q scores [https://doi.org/10.1093/bioinformatics/btv401].

