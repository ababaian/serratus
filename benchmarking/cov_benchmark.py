import argparse
import csv
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='Generates benchmarking statistics for divergence tests.\nWrites TP/FN/FP/TN values in output file.')
parser.add_argument('input',
                    help='FASTA file containing cov genomes.')
parser.add_argument('output',
                    help='Destination filepath for summary CSV file.')
parser.add_argument('--prop_pos', type=float, default=0.05,
                    help='Proportion of genome to mutate for positive control sequence.\nDefault 0.05 (5%% divergence).')
parser.add_argument('--prop_neg', type=float, default=1,
                    help='Proportion of genome to mutate for negative control sequence.\nDefault 1 (100%% divergence).')
args = parser.parse_args()

prop_pos = args.prop_pos
prop_neg = args.prop_neg
cov_genomes_fa = args.input
summary_file = args.output

# create temp dir
tmp_dir = 'tmp_fasta_index'
os.mkdir(tmp_dir)


def write_fastas(record_fwd, fa_fwd, fa_rev):
    """Write forward and reverse non-complement FASTA files."""
    with open(fa_fwd, 'w') as f:
        f.write(record_fwd.format('fasta'))

    seq = str(record_fwd.seq)
    rev_name = f'REVERSE_{record_fwd.id}'
    record_rev = SeqRecord(Seq(seq[::-1]),
                           id=rev_name,
                           name=rev_name,
                           description=rev_name)
    with open(fa_rev, 'w') as f:
        f.write(record_rev.format('fasta'))


def call_msbar(count, sequence, outseq):
    """Call msbar to create divergent sequences. Parameters pre-set."""
    subprocess.check_output(['msbar',
                             '-point', '4',
                             '-block', '0',
                             '-codon', '0',
                             '-count', str(count),
                             '-sequence', sequence,
                             '-outseq', outseq], stderr=subprocess.DEVNULL)


def call_art_illumina(in_seq, out_seq):
    """Call art_illumina to simulate reads. Parameters pre-set."""
    subprocess.check_output(['art_illumina',
                             '--seqSys', 'HS20',
                             '--paired',
                             '--in', in_seq,
                             '--len', '100',
                             '--mflen', '300',
                             '--sdev', '1',
                             '--fcov', '50',
                             '--rndSeed', '666',
                             '--out', out_seq,
                             '--noALN'], stderr=subprocess.DEVNULL)


def get_alignments(target_index, fq_sim_prefix):
    """Call bowtie2 for alignment. Return set of aligned reads."""
    call_bt2 = subprocess.Popen(['bowtie2',
                                 '--very-sensitive-local',
                                 '-x', target_index,
                                 '-1', f'{fq_sim_prefix}1.fq',
                                 '-2', f'{fq_sim_prefix}2.fq'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)
    call_samtools = subprocess.Popen(['samtools', 'view', '-G', '12', '-'],
                                     stdin=call_bt2.stdout, stdout=subprocess.PIPE)
    out_awk = subprocess.check_output(['awk', '!seen[$1]++ {print $1}'],
                                      stdin=call_samtools.stdout)
    return set(filter(None, out_awk.decode().split('\n')))


def get_alignment_stats(record):
    # write fasta and index files

    fa_fwd = os.path.join(tmp_dir, 'fwd.fa')
    fa_rev = os.path.join(tmp_dir, 'rev.fa')
    write_fastas(record, fa_fwd, fa_rev)

    # bowtie2-build genomes

    index_fwd = os.path.join(tmp_dir, 'fwd.index')
    index_rev = os.path.join(tmp_dir, 'rev.index')
    subprocess.check_output(['bowtie2-build', fa_fwd, index_fwd], stderr=subprocess.DEVNULL)
    subprocess.check_output(['bowtie2-build', fa_rev, index_rev], stderr=subprocess.DEVNULL)

    # calculate number of mutations

    genome_len = len(record.seq)
    n_mutations_pos = int(prop_pos * genome_len)
    n_mutations_neg = int(prop_neg * genome_len)

    # simulate divergent genomes

    fa_sim_pos = os.path.join(tmp_dir, 'sim_pos.fa')
    fa_sim_neg = os.path.join(tmp_dir, 'sim_neg.fa')
    call_msbar(n_mutations_pos, fa_fwd, fa_sim_pos)
    call_msbar(n_mutations_neg, fa_fwd, fa_sim_neg)

    # simulate reads from divergent genomes

    fq_sim_pos_prefix = os.path.join(tmp_dir, 'sim_pos_')
    fq_sim_neg_prefix = os.path.join(tmp_dir, 'sim_neg_')
    call_art_illumina(fa_sim_pos, fq_sim_pos_prefix)
    call_art_illumina(fa_sim_neg, fq_sim_neg_prefix)

    # get number of reads (same for pos/neg)
    with open(f'{fq_sim_pos_prefix}1.fq') as fq:
        n_reads = sum(1 for read in FastqGeneralIterator(fq))

    # bowtie2 map, get alignments
    aligned_sim_pos_cov = get_alignments(index_fwd, fq_sim_pos_prefix)
    aligned_sim_neg_cov = get_alignments(index_fwd, fq_sim_neg_prefix)
    aligned_sim_pos_rev_cov = get_alignments(index_rev, fq_sim_pos_prefix)
    aligned_sim_neg_rev_cov = get_alignments(index_rev, fq_sim_neg_prefix)

    tp = len(aligned_sim_pos_cov - aligned_sim_pos_rev_cov)
    fn = n_reads - tp
    fp = len(aligned_sim_neg_cov - aligned_sim_neg_rev_cov)
    tn = n_reads - fp

    return (tp, fn, fp, tn)


if __name__ == '__main__':
    with open(summary_file, 'w') as f:
        output_writer = csv.writer(f, delimiter=',')
        output_writer.writerow(['GENOME', 'TP', 'FN', 'FP', 'TN'])

        for record in SeqIO.parse(cov_genomes_fa, "fasta"):
            # skip REVERSEs
            if record.id.startswith('REVERSE'):
                break

            print(f'Simulating for {record.id}...')
            tp, fn, fp, tn = get_alignment_stats(record)
            output_writer.writerow([record.id, tp, fn, fp, tn])

    # TODO: remove tmp_dir
