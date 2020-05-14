# Serratus Scheduler
This is a scheduler written in Flask, which coordinates the actions of Serratus.

## Development

For local development, I recommend running it locally on your laptop, with the development server.  To do that, from the `serratus/containers` directory, run:

    ./build_containers.sh serratus-scheduler
    podman run -d --rm -p5000:5000 -v $(realpath serratus-scheduler/flask_app):/opt/flask_app:z --name scheduler serratus-scheduler flask run -p 5000

This will run the development server using port 5000 on the local system.  `-v` tells podman to map your local development directory on top of the scheduler's files, so the development server will reload whenever you make a change to your local files.  You can check it's working by watching `podman logs -f scheduler` while making a change to any of the Python files under `serratus-scheduler/flask_app`.

## API Calls (normal use)

### Scheduler Configuration

Get the current server config:

    curl localhost:8000/config > serratus-config.json

This will output a configuration JSON file, like the following:

    {
      "ALIGN_SCALING_CONSTANT": 0.001,
      "ALIGN_SCALING_ENABLE": false,
      "ALIGN_SCALING_MAX": 20,
      "CLEAR_INTERVAL": 300,
      "DL_SCALING_CONSTANT": 0.001,
      "DL_SCALING_ENABLE": false,
      "DL_SCALING_MAX": 10,
      "GENOME": "cov2r",
      "MERGE_SCALING_CONSTANT": 0.001,
      "MERGE_SCALING_ENABLE": false,
      "MERGE_SCALING_MAX": 2,
      "SCALING_INTERVAL": 30
    }

Change these values to what you want, and then update the configuration with another `curl` call:

    curl -T serratus-config.json localhost:8000/config

Notes

 * If you delete keys, they will keep their old values.
 * There should be a comma on every line except the last one.
 * Be careful with types.  true/false fields in particular must *not* be quoted.  In Python, non-null
   strings like "false", "no", and "disabled" will all be converted to True.

### Add SraRunInfo.csv

Send a csv file from SRA to the scheduler, this is the main way to tell Serratus to do work:

    curl -T /path/to/sra_run_info.csv -X POST /jobs/add_sra_run_info/

### View HTML Table

This page shows a table view of all the jobs currently known to the scheduler.

    GET /jobs/

## API Calls (development / debugging)

### Trigger scheduler to check the instance list

This trigger uses Ec2DescribeInstances and checks the instance list against all jobs that are currently in a running state.  If the instance is missing (probably due to being terminated or shut down), it resets the job to its prior state.

    PUT /jobs/clear_terminated

This will now happen automatically, and this request should no longer be needed.

### Get metrics

Use this API to view the Prometheus metrics.  This is mainly used by our Grafana server to create beautiful plots.

    GET /metrics

### Start and End jobs

We use POST to start jobs, and PUT to end jobs.  You can also use these APIs to manually override behaviours in the scheduler.  The general idea is to use `POST /jobs/<type>/` to start a job (ID is returned as part of the job and automatically marked as running), and `PUT /jobs/<type>/<id>` to change its state. For example

#### Reset accession download
This will reset the accession 54 to the "new" state.
```
curl -X POST "localhost:8000/jobs/split/54?state=new&N_paired=0&N_unpaired=0"
```
#### Reset a range of align blocks
This will reset blocks 205-230 to the new state.
```
for BLOCK_ID in $(seq 205 230)
do
  curl -X POST -s "localhost:8000/jobs/align/$BLOCK_ID?state=new"
done
```



