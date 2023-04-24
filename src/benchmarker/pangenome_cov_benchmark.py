import argparse
import csv
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='Run cov_benchmark.py for each genome sequence in a pangenome using defaults and --output_proportions.')
parser.add_argument('pangenome',
                    help='FASTA file containing one record per genome sequence.')
parser.add_argument('out',
                    help='Output CSV file.')
args = parser.parse_args()


script_path = "cov_benchmark.py"


def run_benchmarker(params):
    params_list = ["python", script_path]
    params_list.extend(params)
    run_cmd = " ".join(params_list)
    print(f"Running: {run_cmd}")
    return subprocess.check_output(params_list)



with open(args.out, 'w') as csv_file:
    summary_out = csv.writer(csv_file, delimiter=',')
    summary_out.writerow(['GENOME', 'TP', 'FN', 'FP', 'TN'])

    for record in SeqIO.parse(args.pangenome, "fasta"):
        record_fa = f'{record.id}.fa'
        with open(record_fa, 'w') as f:
            f.write(record.format('fasta'))

        params = [record_fa, '--output_proportions']

        output = run_benchmarker(params)
        output = output.decode().rstrip()
        tp, fn, fp, tn = output.split(",")
        summary_out.writerow([record.id, tp, fn, fp, tn])
        os.remove(record_fa)
