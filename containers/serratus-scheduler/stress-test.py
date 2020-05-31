#!/usr/bin/env python
"""Stress tester for the scheduler"""
from multiprocessing import Process, Queue, Manager
import random
import string
import sys
import queue
import time

import requests

BLOCKS = 10
THREADS = 100


def load_sch(sch, csv, lines):
    """Load some test data into the scheduler"""
    # Reset db
    requests.put("http://{}/db".format(sch))
    data = "\n".join(open(csv).readlines()[: lines + 1])
    requests.post("http://{}/jobs/add_sra_run_info/data.csv".format(sch), data=data)


def serratus_worker(sch, q, type_, params=()):
    """Act in place of worker.sh"""
    start_params = {
        "worker_id": "".join(random.choice(string.ascii_lowercase) for _ in range(10)),
    }
    end_params = {"state": "done"}
    end_params.update(params)

    while True:
        # Get a job
        res = requests.post(
            "http://{}/jobs/{}/".format(sch, type_), params=start_params
        ).json()

        if res["action"] == "wait":
            break

        job_id = res["id"]
        q.append(res["acc_id"])

        # Do nothing, and send the response.
        res2 = requests.post(
            "http://{}/jobs/{}/{}".format(sch, type_, job_id), params=end_params
        )


def main():
    # Load the scheduler
    sch, csv, lines = sys.argv[1:]

    load_sch(sch, csv, int(lines))

    manager = Manager()
    dl_jobs = manager.list()
    al_jobs = manager.list()
    mg_jobs = manager.list()

    argses = (
        (sch, dl_jobs, "dl", {"N_unpaired": BLOCKS}),
        (sch, al_jobs, "align"),
        (sch, mg_jobs, "merge"),
    )

    # Create fake instances of all 3 types.
    procs = []
    for args in argses:
        start = time.monotonic()
        for _ in range(THREADS):
            p = Process(target=serratus_worker, args=args)
            p.start()
            procs.append(p)

        for p in procs:
            p.join()
        print("Total {} time: {}s".format(args[2], time.monotonic() - start))

    # Uniqueness and completeness
    assert len(set(dl_jobs)) == int(lines)
    assert len(set(dl_jobs)) == len(dl_jobs)
    assert len(set(mg_jobs)) == len(mg_jobs)

    # Each acc_id should appear N times for align.
    assert len(al_jobs) == BLOCKS * len(dl_jobs)

    # Same jobs in all parts of the pipeline.
    assert set(dl_jobs) == set(al_jobs)
    assert set(dl_jobs) == set(mg_jobs)


if __name__ == "__main__":
    main()
