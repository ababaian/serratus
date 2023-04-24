import subprocess

script_path = "../../src/benchmarker/cov_benchmark.py"


def run_benchmarker(params):
    params_list = ["python", script_path]
    params_list.extend(params)
    run_cmd = " ".join(params_list)
    print(f"Running: {run_cmd}")
    return subprocess.check_output(params_list)


def test_seq_params():
    pass


def test_prop_params():
    pass


def test_read_set_params():
    params = [
        "local/NC_045512v2r.fa",
        "-v",
        "--pos_reads", "local/sim_pos_",
        "--neg_reads", "local/sim_neg_",
    ]
    output = run_benchmarker(params)
    output = output.decode()
    assert output == "6650,825,0,7475\n"


def test_bowtie2_very_sensitive_local():
    params = [
        "local/NC_045512v2r.fa",
        "-v",
        "--pos_reads", "local/sim_pos_",
        "--neg_reads", "local/sim_neg_",
        "--bowtie2_params=--very-sensitive-local",
    ]
    output = run_benchmarker(params)
    output = output.decode()
    assert output == "7450,25,0,7475\n"


def test_reads_src_params():
    params = [
        "local/NC_045512v2r.fa",
        "-v",
        "--pos_reads_src", "local/pos_reads_src.fa",
        "--neg_reads_src", "local/neg_reads_src.fa",
        "--bowtie2_params=--very-sensitive-local",
        '--art_illumina_params="--rndSeed 666"',
    ]
    output = run_benchmarker(params)
    output = output.decode()
    assert output == "7468,7,0,7475\n"


def test_align_seq_params():
    pass


def test_output_proportions_param():
    pass
