# cov_benchmark.py

## Usage

    $ python cov_benchmark.py -h
    usage: cov_benchmark.py [-h] [--prop_pos PROP_POS] [--prop_neg PROP_NEG] input output
    
    Generates benchmarking statistics for divergence tests.
    Writes TP/FN/FP/TN values in output file.
    
    positional arguments:
      input                FASTA file containing cov genomes.
      output               Destination filepath for summary CSV file.
    
    optional arguments:
      -h, --help           show this help message and exit
      --prop_pos PROP_POS  Proportion of genome to mutate for positive control sequence.
                           Default 0.05 (5% divergence).
      --prop_neg PROP_NEG  Proportion of genome to mutate for negative control sequence.
                           Default 1 (100% divergence).

### Example

    python cov_benchmark.py cov1r.fa alignment_stats.csv

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
