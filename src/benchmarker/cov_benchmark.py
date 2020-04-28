import argparse
import os
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='Runs divergence tests and generates summary statistics.\nOutputs values in format: "TP,FN,FP,TN"')
parser.add_argument('input_seq',
                    help='Single-record FASTA containing input sequence. Used for pos/neg read sets and read alignment, unless other parameters specified.')
parser.add_argument('--pos_seq', metavar='FASTA',
                    help='FASTA sequences used for simulation of positive read set.\nDefault: input sequence.')
parser.add_argument('--neg_seq', metavar='FASTA',
                    help='FASTA sequences used for simulation of negative read set.\nDefault: reverse non-complement of input sequence.')
parser.add_argument('--prop_pos', metavar='PROP', type=float,
                    help='Proportion of pos_seq to mutate for simulation of positive read set.\nDefault: 0.05 (5%% divergence).')
parser.add_argument('--prop_neg', metavar='PROP', type=float,
                    help='Proportion of neg_seq to mutate for simulation of negative read set.\nDefault: 0.05 (5%% divergence).')
parser.add_argument('--pos_reads', metavar='FQ_PREFIX',
                    help='Positive read set. Paired-end files are {{FQ_PREFIX}}1.fq and {{FQ_PREFIX}}2.fq.\nIf specified, ignores pos_seq and prop_pos.\nDefault: reads derived from pos_seq with prop_pos applied.')
parser.add_argument('--neg_reads', metavar='FQ_PREFIX',
                    help='Negative read set. Paired-end files are {{FQ_PREFIX}}1.fq and {{FQ_PREFIX}}2.fq.\nIf specified, ignores neg_seq and prop_neg.\nDefault: reads derived from neg_seq with prop_neg applied.')
parser.add_argument('--pos_align_seq', metavar='FASTA',
                    help='Reference sequence for read alignment.\nReads mapped to this SEQ only will be classified as TP/FP.\nDefault: input sequence.')
parser.add_argument('--neg_align_seq', metavar='FASTA',
                    help='Reference sequence for read alignment.\nReads mapped to this SEQ will be classified as TN/FN, along with unmapped reads.\nDefault: reverse non-complement of input sequence.')
parser.add_argument('--msbar_params', metavar='STR',
                    help='Additional parameters to pass to the msbar command.')
parser.add_argument('--art_illumina_params', metavar='STR',
                    help='Additional parameters to pass to the art_illumina command.')
parser.add_argument('--bowtie2_params', metavar='STR',
                    help='Additional parameters to pass to the bowtie2 command.')
parser.add_argument('-v', action='store_true',
                    help='Additional parameters to pass to the bowtie2 command.')
args = parser.parse_args()


# create temp dir
tmp_dir = 'tmp'
os.makedirs(tmp_dir, exist_ok=True)


class Command(object):
    def __init__(self, name, args):
        self.name = name
        self.args = args

    def run(self, quiet=True, pipe=False):
        args_list = [self.name]
        args_list.extend([str(val) for pair in self.args.items()
                          for val in pair if val is not None])
        subprocess.check_output(args_list, stderr=(subprocess.DEVNULL if quiet else None))
        # TODO: implement pipe


def printv(msg):
    if args.v:
        print(msg, file=sys.stderr)


def write_reverse_fasta(record, fa_out):
    """Write reverse non-complement FASTA file."""
    seq = str(record.seq)
    rev_name = f'REVERSE_{record.id}'
    record_rev = SeqRecord(Seq(seq[::-1]),
                           id=rev_name,
                           name=rev_name,
                           description=rev_name)
    with open(fa_out, 'w') as f:
        f.write(record_rev.format('fasta'))


def get_alignments(target_index, fq_sim_prefix):
    """Call bowtie2 for alignment. Return set of aligned reads."""
    call_bt2 = subprocess.Popen(['bowtie2',
                                 '--very-sensitive-local',
                                 '-x', target_index,
                                 '-1', f'{fq_sim_prefix}1.fq',
                                 '-2', f'{fq_sim_prefix}2.fq'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)
    call_samtools = subprocess.Popen(['samtools', 'view', '-G', '12', '-'], stdin=call_bt2.stdout, stdout=subprocess.PIPE)
    out_awk = subprocess.check_output(['awk', '!seen[$1]++ {print $1}'], stdin=call_samtools.stdout)
    return set(filter(None, out_awk.decode().split('\n')))


def simulate_read_set(seq, prop, reads_prefix):
    """Simulate read set given sequence, mutation proportion, and prefix for output reads."""
    seq_len = sum(len(r.seq) for r in SeqIO.parse(seq, "fasta"))
    sim_fa = os.path.join(tmp_dir, 'sim.fa')
    n_mutations = int(prop * seq_len)
    args_msbar = {
        '-point': 4,
        '-block': 0,
        '-codon': 0,
        '-count': n_mutations,
        '-sequence': seq,
        '-outseq': sim_fa
    }
    cmd_msbar = Command('msbar', args_msbar)
    cmd_msbar.run()
    args_art_illumina = {
        '--seqSys': 'HS20',
        '--paired': None,
        '--in': sim_fa,
        '--len': 100,
        '--mflen': 300,
        '--sdev': 1,
        '--fcov': 50,
        '--rndSeed': 666,
        '--out': reads_prefix,
        '--noALN': None
    }
    cmd_art_illumina = Command('art_illumina', args_art_illumina)
    cmd_art_illumina.run()
    os.remove(sim_fa)


def get_alignment_stats(tmp_dir):

    # get input SeqRecord
    record = next(SeqIO.parse(args.input_seq, "fasta"))

    reverse_written = False

    if not args.pos_reads:
        printv('Simulating positive read set...')

        if not args.pos_seq:
            printv('pos_seq not specified - using input_seq.')
            args.pos_seq = args.input_seq
        if not args.prop_pos:
            printv('prop_pos not specified - using 0.05.')
            args.prop_pos = 0.05

        args.pos_reads = os.path.join(tmp_dir, 'sim_pos_')
        simulate_read_set(args.pos_seq, args.prop_pos, args.pos_reads)
    else:
        printv(f'Using positive read set files {args.pos_reads}(1,2).fq.')
        if args.pos_seq:
            printv(f'WARNING: pos_reads already specified - pos_seq={args.pos_seq} will be ignored.')
        if args.prop_pos:
            printv(f'WARNING: pos_reads already specified - prop_pos={args.prop_pos} will be ignored.')


    if not args.neg_reads:
        printv('Simulating negative read set...')

        if not args.neg_seq:
            printv('neg_seq not specified - using reverse non-complement of input_seq.')
            # write reverse FASTA
            args.neg_seq = os.path.join(tmp_dir, 'neg_seq.fa')
            write_reverse_fasta(record, args.neg_seq)
            reverse_written = True
        if not args.prop_neg:
            printv('prop_neg not specified - using 0.05.')
            args.prop_neg = 0.05

        args.neg_reads = os.path.join(tmp_dir, 'sim_neg_')
        simulate_read_set(args.neg_seq, args.prop_neg, args.neg_reads)
    else:
        printv(f'Using negative read set files {args.neg_reads}(1,2).fq.')
        if args.neg_seq:
            printv(f'WARNING: neg_reads already specified - neg_seq={args.neg_seq} will be ignored.')
        if args.prop_neg:
            printv(f'WARNING: neg_reads already specified - prop_neg={args.prop_neg} will be ignored.')


    if not args.pos_align_seq:
        printv('pos_align_seq not specified - using input_seq.')
        args.pos_align_seq = args.input_seq
    else:
        printv(f'Using pos_align_seq={args.pos_align_seq}.')
    if not args.neg_align_seq:
        printv('neg_align_seq not specified - using reverse non-complement of input_seq.')
        if reverse_written:
            args.neg_align_seq = args.neg_seq
        else:
            args.neg_align_seq = os.path.join(tmp_dir, 'neg_align_seq.fa')
            write_reverse_fasta(record, args.neg_align_seq)
    else:
        printv(f'Using neg_align_seq={args.neg_align_seq}.')

    printv('Aligning read sets...')

    # build FASTA indexes
    index_pos = os.path.join(tmp_dir, 'pos.index')
    subprocess.check_output(['bowtie2-build', args.pos_align_seq, index_pos], stderr=subprocess.DEVNULL)
    index_neg = os.path.join(tmp_dir, 'neg.index')
    subprocess.check_output(['bowtie2-build', args.neg_align_seq, index_neg], stderr=subprocess.DEVNULL)

    # calculate statistics

    with open(f'{args.pos_reads}1.fq') as fq:
        n_reads_pos = sum(1 for r in FastqGeneralIterator(fq))
    with open(f'{args.neg_reads}1.fq') as fq:
        n_reads_neg = sum(1 for r in FastqGeneralIterator(fq))

    aligned_sim_pos_cov = get_alignments(index_pos, args.pos_reads)
    aligned_sim_neg_cov = get_alignments(index_pos, args.neg_reads)
    aligned_sim_pos_rev_cov = get_alignments(index_neg, args.pos_reads)
    aligned_sim_neg_rev_cov = get_alignments(index_neg, args.neg_reads)

    tp = len(aligned_sim_pos_cov - aligned_sim_pos_rev_cov)
    fn = n_reads_pos - tp
    fp = len(aligned_sim_neg_cov - aligned_sim_neg_rev_cov)
    tn = n_reads_neg - fp

    return (tp, fn, fp, tn)


if __name__ == '__main__':
    tp, fn, fp, tn = get_alignment_stats(tmp_dir)
    print(f'{tp},{fn},{fp},{tn}')
