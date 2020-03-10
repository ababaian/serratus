#!/usr/bin/bash
# Curl the server 1000 times, with 10 threads, to check it doesn't produce
# invalid data.
N=100
THREADS=10

function go {
    for i in $(seq 1 $N); do
        curl -s http://localhost:8000/start_split_job | grep acc_id
    done
}

for i in $(seq 1 "$THREADS"); do
    rm -f out."$i".log
    go >> out."$i".log &
done
wait
