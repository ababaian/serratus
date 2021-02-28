# Gunicorn and Postgres Tuning

We're switching Gunicorn to async workers, which has some implications across
the application.  It's possible that our application will be CPU-bound (parsing
JSON), or IO-bound (waiting for Postgres), so let's make sure it behaves
appropriately in both cases.

## The app

I wrote a test app with two routes.  `/io` which is clearly IO-bound:

    @app.route("/io")
    def io():
        with db.engine.connect() as con:
            con.execute("SELECT pg_sleep(5);")
        return "hello there!\n"

and `/cpu`, which is clearly CPU-bound:

    @app.route("/cpu")
    def cpu():
         return str(sum(range(100000000))) # About 2 seconds on my laptop.

To test configurations, I run the server using Gunicorn, and spin up 3000 curls.
I'm intentionally overloading the server, so I expect the results to trickle in.
Unacceptable are server crashes, invalid responses, and so on.

## SQLAlchemy QueuePools

Let's run gunicorn with the following settings:

    worker_class = "gevent"
    workers = 1

Right away, when we run our "1000 curls", we see problems:

    $ for i in $(seq 1 3000); do curl localhost:8001/io & done; time wait
    QueuePool limit of size 5 overflow 10 reached, connection timed out, timeout 30.

This is a problem, and it gives us some hints about how SQLAlchemy manages connections, explained very well by [http://www.javaear.com/question/57844921.html], and in more detail in the documentation: https://duckduckgo.com/?t=ffab&q=sqlalchemy+connection+pooling&atb=v143-1&ia=web

To summarize, SQLAlchemy maintains a pool of open connections.  This connection pool starts at zero connections, and grows until it reaches `pool_size` (5 by default).  If all 5 connections are occupied, it creates up to `overflow` extra connections, but closes them as soon as they're returned.

## CPU Load

Let's run the same curl command as above, but with the `/cpu` endpoint, to see how this would handle being CPU-bound.  In this case we find that all of our CPUs are initially maxxed out, but the processes drop out, and in the end, only one process is handling the requests.  This is also a problem.

## Gunicorn Worker Connections

It turns out these have a common cause: gunicorn by default allows each worker to process up to 1000 requests.  https://docs.gunicorn.org/en/stable/settings.html#worker-connections

This causes problems with CPU-bound loads, because a worker can pick up to 1000 jobs, but then if those jobs are CPU-heavy, it doesn't actually have the capacity to process them all in a single process.  Some of them will eventually timeout, creating stability problems.

This is the same situation for IO-heavy loads.  The 15-connection limit puts a limit on the amount of throughput a process can do.  If it tries to handle 1000 requests through 15 connections, those requests can time out and fail.

Let's try reducing this to a more sane value.  It should be somewhere close to the maximum number of open connections.  If it is lower, we'll be wasting connections.  If it's higher, our requests have the potential to block waiting for a connection.  In our case, we have a mix of high latency operations (merge queries) and low latency operations (align queries), which differ by around two orders of magnitude in terms response time.  I don't want the fast operations to get stuck in a queue with the slow operations, so I've set `worker_connections` less than `pool_size + overflow`.

## Postgres Max Connections

The last thing to tune is postgres `max_connections`.  This is fairly easy to calculate, by `worker_connections * workers`.  For a 32-CPU m5.8xlarge instance, we'll have `33 * 15 = 495` connections maximum from Flask.  Let's set this to 1000, so that we have some overhead room.

I've found that just increasing `max_connections` in Postgres (version 12) works without any other steps.  I didn't need to change kernel parameters or anything complicated.

