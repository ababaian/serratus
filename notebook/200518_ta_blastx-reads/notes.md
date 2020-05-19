# Introduction

Based on a hypothesis of Robert's, I ran the reads of sample
ERR2756788 (e.g., "Frankie") against the proteome of all of the
sequences hit by Bowtie2. Robert will analyze this output to see if we
can come up with ultra-sensitive contig assembly or selectin
(depending on which kind of assembly one does.


## Methods

Download the ERR2756788 summary file, scrape out the NCBI accessions,
use EDirect to convert into proteins, and download the protein
sequences:
```
time make -f notebook/200518_ta_tblastx-reads.make /media/storage/tblastx/cov-hit-prots.fsa
```

Generate a BLAST v5 DB from the multi-FASTA protein file:
```
time make -f notebook/200518_ta_tblastx-reads.make /media/storage/tblastx/cov-hit-prots.fsa.pdb
```

Convert the FASTQ reads into FASTA format, and concatenate the paired
end files into a single multi-FASTA file:

```
date; time make -f notebook/200518_ta_tblastx-reads.make /media/storage/tblastx/query-reads.fa; echo $?; date
```

Run blastx in a massively-parallel fashion:

```
date; time make -f notebook/200518_ta_tblastx-reads.make parallel-blastx; echo $?; date
```

Stage the data on S3:

```
date; time make -f notebook/200518_ta_tblastx-reads.make stage-on-s3;
echo $?; date
```

