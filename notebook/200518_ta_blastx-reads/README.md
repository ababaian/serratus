# Introduction

This folder contains various experiments related to how to search
prospective SRA runs quickly to determine whether there is CoV-related
reads which warrant more detailed analysis.

## SRA BlastN_

A poorly-documented feature of the `sra-toolkit` is a program called
`blastn_vdb`. It is a binary that allows the NCBI Blast suite
subroutines to access sequence data in SRA-formatted files directly,
without a need to first use a tool like `fastq-dump` to export the
data to a flat-file format. This provides immediate savings in
processing time if the search itself doesn't eat up those savings.

### Methodology

I sought to evaluate the performance of using `blastn_vdb` both in
terms of resource utilization (i.e., CPUs, RAM, etc.) and
sensitivity. I repurposed some simulated reads that Artem had used in
the past to evaluate the fitness of using Bowtie2 for detecting
reads from divergent Coronavirus genomes:

`https://github.com/ababaian/serratus/blob/master/notebook/200411_CoV_Divergence_Simulations.ipynb`

#### Testing Sensitivity

I wanted to reuse Artem's testing set-up, so I had to use `blastn`
from the NCBI Blast package instead of testing `blastn_vdb` directly,
as it can only work on SRA files. The code that does the hit search
and alignment should be exactly the same, so it is used here as a
proxy.

Download & set up some required software:

```shell
mkdir -p third-party
cd third-party
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz
tar xzf ncbi-blast-2.10.0+-x64-linux.tar.gz
export PATH=$PATH:$PWD/third-party/ncbi-blast-2.10.0+/bin
```

First, we fetch the benchmark simulated reads:
```shell
mkdir -p benchmark
time aws s3 sync s3://serratus-public/notebook/200411/fq fq
cd fq
for file in `ls *.gz`; do gunzip $file; done
```

Next, we build Blast DBs for each one:
```shell
cd benchmark/fq
MU=(0 30 300 1500 3000 4500 6000 7500 9000 10500 12000)
for mu in ${MU[@]}
do
	mkdir -p sim.cov.${mu}
	pushd sim.cov.${mu}
	cat <(sed -n '1~4s/^@/>/p;2~4p' ../sim.cov.${mu}_1.fq) \
	    <(sed -n '1~4s/^@/>/p;2~4p' ../sim.cov.${mu}_2.fq) \
	    > sim.cov.${mu}.fasta
	makeblastdb -dbtype nucl -in sim.cov.${mu}.fasta
	popd
done	
```

Finally, we run `blastn` with the SARS-CoV-2 reference genome
`NC_045512.2` as the query, and the indexed reads as the subject
(database), mirroring the orientation used by `blastn_vdb`:

```shell
cd benchmark
mkdir -p blastn-output
## 
MU=(0 30 300 1500 3000 4500 6000 7500 9000 10500 12000)
rm -f blastn-output/alignment-rate.txt
touch blastn-output/alignment-rate.txt
for mu in ${MU[@]}
do
        time blastn \
            -db fq/sim.cov.${mu}/sim.cov.${mu}.fasta \
            -query ../SARS-CoV-2.fa \
            -out blastn-output/sim.cov.${mu}.tsv \
            -evalue 0.001 \
            -max_target_seqs 100000 \
            -task blastn \
            -outfmt "7 std qlen slen"
    awk -v tot="`egrep -c "^>" fq/sim.cov.${mu}/sim.cov.${mu}.fasta`" \
    	'!/^#/ && !a[$2]++ { count++}END{print count/tot}' \
	blastn-output/sim.cov.${mu}.tsv \
	>> blastn-output/alignment-rate.txt
done
```

#### Testing Performance

I wrote a Snakemake file to automate the testing of `blastn_vdb` (see
`Snakefile` in this directory). I tested its performance against the
Frankie SRA run ('ERR2756788'), a SARS-CoV-2 positive (CoV+) SRA run
('SRR11454614'), and a SARS-CoV-2 negative (CoV-) SRA run
('ERR3568641'), inspired by Robert's benchmarking set-up as found
here:

`https://github.com/ababaian/serratus/blob/54856f0f86ff8f0af3ed6f9e87f0d7d04819c570/notebook/2000605_rce_diamond.pdf`

### Results

#### Sensitivity

Here is a table of the alignment rate of `blastn` as a function of the
mutation rate:


Mutation % | Align %
--- | ---
0 | 1
0.001 | 1
0.01 | 1
0.05 | 0.999866
0.1 |  0.999799
0.15 | 0.994983
0.2 |  0.980936
0.25 | 0.946154
0.3 |  0.856455
0.35 | 0.753612
0.4 |  0.632642


Here is a plot of the alignment rate of `blastn` as a function of the
mutation rate:

![Sensitivity plot](ali)

#### Performance

Running `blastn_vdb` on a pre-fetched SRA file took only one minute
using a single thread, using the SARS-CoV-2 genome as the
"query". The runtime seemed to increase linearly with the size of the
query, so that putting in another SARS-CoV-2 genome would double the
run-time. Memory consumption was negligible.


### Discussion & Future Work

The use of `blastn_vdb` would allow for a significant improvement in
the Serratus pipeline, as Artem reports that there's roughly a 1:10
ratio between the time it takes to pre-fetch an SRA file, and the time
taken to dump out the FASTQ sequence to disk. The sensitivity is
superior to that of Bowtie2: while Bowtie2's alignment rate bottoms
out to near zero at 40% divergence, BlastN is still aligning at 63%.

The only question is about `blastn_vdb`'s ability to handle a larger
reference database as its "query", as currently it wouldn't be
performant taking the entire `cov3ma` database as-is. One of the
drivers for making such an expansive cov3ma database was the rationale
that if Bowtie2's sensitivity dropped off, that we can compensate by
having a greater diversity in the reference database. With BlastN's
superior sensitivity, it might be a wash. Further testing would
determine where the break-even point might be. One thing that would
bear immediate fruit would be to make the cov3ma database
non-redundant, as currently has many copies of complete coronavirus
genomes. An approach like what is taken by `Centrifuge` to remove
redundancy in a genome database would be useful for both Bowtie2 or
BlastN:

`http://genome.cshlp.org/content/early/2016/11/16/gr.210641.116`

Another use of `blastn_vdb` would be as a "predicate", to quickly
determine in one minute whether a given sample is worth dumping out as
FASTQ to disk, and checking more carefully with a larger reference
database.

Due to a limitation in `blastn_vdb`, it can only run
single-threaded. This limitation can be circumvented by using GNU
Parallel to split up the reference database of CoV sequences (which,
according to `blastn_vdb`, are the query sequences) , running a single
sequence from the CoV sequences against the SRA file at a time.

Currently only the `blastn` program is available in the SRA toolkit
for aligning against SRA files directly. We can simulate a
`blastx`-like program by, for a specificed translation code, render
each protein sequence in all of the possible translations from peptide
sequence back to nucleotide sequence. Using informative identifiers,
these back-translated sequences can be mapped to the original protein
sequences using a down-stream script. Since we are using BlastN, we
can use the IUPAC ambiguous nucleotide codes to minimize the number of
nucleotide sequences created for a given peptide sequence.

When run on a positive control SRA run (SRR11454614; patient with
SARS-CoV-2), it can take a very long time to complete due to the large
number of alignments. A work-around here is to make use of the
`-max_target_seqs` option, which will "short circuit" the search as
soon as the specified number of hits have been found. When run against
the "Frankie" SRA (ERR2756788) with a max target limit of 10,000,000
hits, ~1,800 hits were found. A reasonable upper-limit for the max
number of hits that ensures that all divergent genomes are found, but
limits the runtime for corner-cases where the SRA is a sample with
high-abundance of a known coronavirus, an upper-limit of 100,000 is
reasonable, as it doesn't seem to impact the run-time.

It is unclear whether the "database" (i.e., the one or more SRA files)
are indexed before `blastn_vdb` begins its search. I have put in a
request to the authors for more information.

### Conclusion

Use of `blastn_vdb` to quickly screen SRA files without having to
invest (in terms of time and disk space) in dumping its contents out
to FASTQ files would save significant resources and allow Serratus to
scale to higher throughput. It has greater sensitivity than Bowtie2.
