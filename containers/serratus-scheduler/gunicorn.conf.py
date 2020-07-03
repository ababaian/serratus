import multiprocessing
from prometheus_client import multiprocess as prom_mp

workers = multiprocessing.cpu_count() * 2 + 1


def child_exit(server, worker):
    prom_mp.mark_process_dead(worker.pid)
