# cov_benchmark.py

## Definitions

- `input_seq` : Single-record FASTA containing input sequence. Used for pos/neg read sets and read alignment, unless other parameters specified.
- `pos_seq` : Positive source sequence
- `neg_seq` : Negative source sequence
- `prop_pos` : Proportion of positive source sequence to mutate
- `prop_neg` : Proportion of negative source sequence to mutate
- `pos_reads_src` : Sequence used for simulation of positive reads
- `neg_reads_src` : Sequence used for simulation of negative reads
- `pos_reads` : Positive reads
- `neg_reads` : Negative read
- `pos_align_seq` : Positive alignment reference sequence
- `neg_align_seq` : Negative alignment reference sequence

## Procedure

1. Simulate positive/negative read sets
1. Re-align read sets to positive/negative alignment reference sequences
1. Calculate counts:
    - `TP` : number of reads in `pos_reads` that align to `pos_align_seq` only
    - `FN` : number of reads in `pos_reads` that align to `neg_align_seq` or remain unmapped
    - `FP` : number of reads in `neg_reads` that align to `pos_align_seq` only
    - `TN` : number of reads in `neg_reads` that align to `neg_align_seq` or remain unmapped
1. Convert counts to proportions:
    - `TP` : `TP/n_pos_reads`
    - `FN` : `FN/n_pos_reads`
    - `FP` : `FP/n_neg_reads`
    - `TN` : `TN/n_neg_reads`


## Usage

```
$ python cov_benchmark.py -h
usage: cov_benchmark.py [-h] [--pos_seq FASTA] [--neg_seq FASTA] [--prop_pos PROP] [--prop_neg PROP] [--pos_reads_src FASTA] [--neg_reads_src FASTA] [--pos_reads FQ_PREFIX]
                        [--neg_reads FQ_PREFIX] [--pos_align_seq FASTA] [--neg_align_seq FASTA] [--msbar_params STR] [--art_illumina_params STR] [--bowtie2_params STR]
                        [--output_proportions] [-v]
                        input_seq

Runs divergence tests and generates summary statistics.
Outputs values in format: "TP,FN,FP,TN"

positional arguments:
  input_seq             Single-record FASTA containing input sequence. Used for pos/neg read sets and read alignment, unless other parameters specified.

optional arguments:
  -h, --help            show this help message and exit
  --pos_seq FASTA       FASTA sequences used for divergence and simulation of positive read set.
                        Default: input sequence.
  --neg_seq FASTA       FASTA sequences used for divergence and simulation of negative read set.
                        Default: reverse non-complement of input sequence.
  --prop_pos PROP       Proportion of pos_seq to mutate for simulation of positive read set.
                        Default: 0.05 (5% divergence).
  --prop_neg PROP       Proportion of neg_seq to mutate for simulation of negative read set.
                        Default: 0.05 (5% divergence).
  --pos_reads_src FASTA
                        FASTA sequences used for simulation of positive read set.
                        Default: pos_seq mutated with divergence defined by prop_pos.
  --neg_reads_src FASTA
                        FASTA sequences used for simulation of negative read set.
                        Default: neg_seq mutated with divergence defined by prop_neg.
  --pos_reads FQ_PREFIX
                        Positive read set. Paired-end files are {{FQ_PREFIX}}1.fq and {{FQ_PREFIX}}2.fq.
                        If specified, ignores pos_seq and prop_pos.
                        Default: reads derived from pos_seq with prop_pos applied.
  --neg_reads FQ_PREFIX
                        Negative read set. Paired-end files are {{FQ_PREFIX}}1.fq and {{FQ_PREFIX}}2.fq.
                        If specified, ignores neg_seq and prop_neg.
                        Default: reads derived from neg_seq with prop_neg applied.
  --pos_align_seq FASTA
                        Reference sequence for read alignment.
                        Reads mapped to this SEQ only will be classified as TP/FP.
                        Default: input sequence.
  --neg_align_seq FASTA
                        Reference sequence for read alignment.
                        Reads mapped to this SEQ will be classified as TN/FN, along with unmapped reads.
                        Default: reverse non-complement of input sequence.
  --msbar_params STR    Additional parameters to pass to the msbar command.
  --art_illumina_params STR
                        Additional parameters to pass to the art_illumina command.
  --bowtie2_params STR  Additional parameters to pass to the bowtie2 command.
  --output_proportions  Output proportions relative to number of reads in set, opposed to counts
  -v                    Enable verbose logging.
```

### Examples

```
$ python cov_benchmark.py NC_045512v2r.fa -v
Simulating positive read set...
pos_reads_src not specified - using NC_045512v2r.fa with prop=0.05 divergence.
Running: msbar -point 4 -block 0 -codon 0 -count 1495 -sequence NC_045512v2r.fa -outseq tmp/pos_reads_src.fa
Running: art_illumina --in tmp/pos_reads_src.fa --out tmp/sim_pos_ --paired --seqSys HS20 --len 100 --mflen 300 --sdev 1 --fcov 50 --noALN
Simulating negative read set...
neg_reads_src not specified - using tmp/neg_seq.fa with prop=0.05 divergence.
Running: msbar -point 4 -block 0 -codon 0 -count 1495 -sequence tmp/neg_seq.fa -outseq tmp/neg_reads_src.fa
Running: art_illumina --in tmp/neg_reads_src.fa --out tmp/sim_neg_ --paired --seqSys HS20 --len 100 --mflen 300 --sdev 1 --fcov 50 --noALN
pos_align_seq not specified - using NC_045512v2r.fa.
neg_align_seq not specified - using reverse non-complement of NC_045512v2r.fa.
Aligning read sets...
Running: bowtie2 -x tmp/pos.index -1 tmp/sim_pos_1.fq -2 tmp/sim_pos_2.fq
Running: bowtie2 -x tmp/pos.index -1 tmp/sim_neg_1.fq -2 tmp/sim_neg_2.fq
Running: bowtie2 -x tmp/neg.index -1 tmp/sim_pos_1.fq -2 tmp/sim_pos_2.fq
Running: bowtie2 -x tmp/neg.index -1 tmp/sim_neg_1.fq -2 tmp/sim_neg_2.fq
6807,668,0,7475
```

Custom mutation proportions for generation of read sequence sources:

```
$ python cov_benchmark.py NC_045512v2r.fa -v \
    --prop_pos 0.02 \
    --prop_neg 0.10
```

Custom read sequence sources:

```
$ python cov_benchmark.py NC_045512v2r.fa -v \
      --pos_reads_src pos_reads_src.fa \
      --neg_reads_src neg_reads_src.fa
```

Custom read sets:

```
$ python cov_benchmark.py NC_045512v2r.fa -v \
      --pos_reads sim_pos_ \
      --neg_reads sim_neg_
```

Custom alignment sequences:

```
$ python cov_benchmark.py NC_045512v2r.fa -v \
      --pos_reads sim_pos_ \
      --neg_reads sim_neg_
```

Custom command parameters:

```
$ python cov_benchmark.py NC_045512v2r.fa -v \
      --art_illumina_params="--rndSeed 666" \
      --bowtie2_params="--very-sensitive-local"
```

## Dependencies

- `python>=3.6`

### Python Packages

- `biopython`

### Command-line Tools

- `bowtie2`
- `msbar`
- `art_illumina`
- `samtools`
- `awk`

## Testing

```
cd test/benchmarker
pytest -s
```
