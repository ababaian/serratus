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


def parse_flag_param_string(paramstring):
    """Return parameter dict dict from argument string.
    Assume each argument key starts with '-' and has one or no values following it.
    Does not work with values that have whitespace within."""
    items = paramstring.split(' ')
    args_dict = {}
    i = 0
    while i < len(items):
        if items[i].startswith('-'):
            if i != len(items) - 1 and not items[i + 1].startswith('-'):
                args_dict[items[i]] = items[i + 1]
                i += 2
            else:
                args_dict[items[i]] = None
                i += 1
    return args_dict


class Command(object):
    def __init__(self, name, positional_params=None, flag_params=None):
        self.name = name
        self.positional_params = positional_params
        self.flag_params = flag_params

    def run(self, quiet=True, pipe_in=None, pipe_out=False):
        params_list = [self.name]
        if self.positional_params:
            params_list.extend([str(val) for val in self.positional_params])
        if self.flag_params:
            params_list.extend([str(val) for pair in self.flag_params.items()
                              for val in pair if val is not None])
        if pipe_out:
            return subprocess.Popen(params_list,
                                    stdin=pipe_in,
                                    stdout=subprocess.PIPE,
                                    stderr=(subprocess.DEVNULL if quiet else None))
        else:
            return subprocess.check_output(params_list,
                                           stdin=pipe_in,
                                           stderr=(subprocess.DEVNULL if quiet else None))

    def update_flag_params(self, paramstring):
        if not paramstring:
            return
        new_params = parse_flag_param_string(paramstring)
        self.flag_params.update(new_params)


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
    bt2_params_default = {
        '-x': target_index,
        '-1': f'{fq_sim_prefix}1.fq',
        '-2': f'{fq_sim_prefix}2.fq'
    }
    cmd_bt2 = Command('bowtie2', flag_params=bt2_params_default)
    cmd_bt2.update_flag_params(args.bowtie2_params)
    call_bt2 = cmd_bt2.run(pipe_out=True)

    samtools_params = {
        '-G': 12,
        '-': None
    }
    cmd_samtools = Command('samtools', positional_params=['view'], flag_params=samtools_params)
    call_samtools = cmd_samtools.run(pipe_in=call_bt2.stdout, pipe_out=True)

    cmd_awk = Command('awk', positional_params=['!seen[$1]++ {print $1}'])
    out_awk = cmd_awk.run(pipe_in=call_samtools.stdout)

    return set(filter(None, out_awk.decode().split('\n')))


def simulate_read_set(seq, prop, reads_prefix):
    """Simulate read set given sequence, mutation proportion, and prefix for output reads."""
    seq_len = sum(len(r.seq) for r in SeqIO.parse(seq, "fasta"))
    sim_fa = os.path.join(tmp_dir, 'sim.fa')
    n_mutations = int(prop * seq_len)
    msbar_params_default = {
        '-point': 4,
        '-block': 0,
        '-codon': 0,
        '-count': n_mutations,
        '-sequence': seq,
        '-outseq': sim_fa
    }
    cmd_msbar = Command('msbar', flag_params=msbar_params_default)
    cmd_msbar.update_flag_params(args.msbar_params)
    cmd_msbar.run()
    art_illumina_params_default = {
        '--in': sim_fa,
        '--out': reads_prefix,
        '--seqSys': 'HS20',
        '--paired': None,
        '--len': 100,
        '--mflen': 300,
        '--sdev': 1,
        '--fcov': 50,
        '--rndSeed': 666,
        '--noALN': None
    }
    cmd_art_illumina = Command('art_illumina', flag_params=art_illumina_params_default)
    cmd_art_illumina.update_flag_params(args.art_illumina_params)
    cmd_art_illumina.run()
    os.remove(sim_fa)


def get_alignment_stats():

    # get input SeqRecord
    record_iter = SeqIO.parse(args.input_seq, "fasta")
    record = next(record_iter)
    if next(record_iter, None) is not None:
        raise ValueError('Input sequence has more than one record. Exiting.')

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
    tp, fn, fp, tn = get_alignment_stats()
    print(f'{tp},{fn},{fp},{tn}')
