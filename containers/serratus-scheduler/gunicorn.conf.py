import multiprocessing
from prometheus_client import multiprocess as prom_mp

# Use async workers, since we spend a lot of time waiting for postgres.
worker_class = "gevent"

# With async workers, we can saturate the CPUs even when some threads are in
# IO.  Use cpu_count() + 1 assuming there's *some* time switching threads or
# getting data.
workers = multiprocessing.cpu_count() + 1

# The default of 1000 is *way* too high.  Just by probability some workers
# will end up processing slower requests, and "fast" requests will get stuck
# waiting in their SQLAlchemy QueuePools behind slow requests.  Keeping the
# number of requests actively being processed by each worker *less than* the
# worker's maximum number of active requests, means that jobs will go to idle
# workers.  The maximum number of requests in progress is defined by SQLAlchemy
# QueuePool.pool_size + QueuePool.overflow (defaults 5, and 10, for a total of
# 15 connections at a time).  This allows for the following mix under heavy load:
#  * 1 in CPU work (eg: parsing JSON)
#  * 9 waiting for Postgres
worker_connections = 15

# Due to the same, we'll need to increase the number of Postgres max_connections,
# up to:
#     workers * (pool_size + overflow)
# For an m5.8xlarge, with 32vCPUs, that's (32 + 1) * (5 + 10) = 500
# We'll use 1000 connections, for good measure.


def child_exit(server, worker):
    prom_mp.mark_process_dead(worker.pid)
