# Serratus Scheduler
This is a scheduler written in Flask, which coordinates the actions of Serratus.

## Development

For local development, I recommend running it locally on your laptop, with the development server.  To do that, from the `serratus/containers` directory, run:

    ./build_containers.sh serratus-scheduler
    podman run -d --rm -p5000:5000 -v $(realpath serratus-scheduler/flask_app):/opt/flask_app:z --name scheduler serratus-scheduler flask run -p 5000

This will run the development server using port 5000 on the local system.  `-v` tells podman to map your local development directory on top of the scheduler's files, so the development server will reload whenever you make a change to your local files.  You can check it's working by watching `podman logs -f scheduler` while making a change to any of the Python files under `serratus-scheduler/flask_app`.

## API Calls

### Scheduler Configuration

Get the current server config:

    GET /config

Update configuratation values, from JSON input.

    PUT /config

For example `curl -T <(echo '{"GENOME": "cov1"}') localhost:8000/config` will update the genome name.

### Add SraRunInfo.csv

Send a csv file from SRA to the scheduler, this is the main way to tell Serratus to do work:

    curl -T /path/to/sra_run_info.csv -X POST /jobs/add_sra_run_info/

### View HTML Table

This page shows a table view of all the jobs currently known to the scheduler.

    GET /jobs/

### Trigger scheduler to check the instance list

This trigger uses Ec2DescribeInstances and checks the instance list against all jobs that are currently in a running state.  If the instance is missing (probably due to being terminated or shut down), it resets the job to its prior state.

    PUT /jobs/clear_terminated

### Get metrics

Use this API to view the Prometheus metrics.  This is mainly used by our Grafana server to create beautiful plots.

    GET /metrics

### Start and End jobs

We use POST to start jobs, and PUT to end jobs.  You can also use these APIs to manually override behaviours in the scheduler.  The general idea is to use `POST /jobs/<type>/` to start a job (ID is returned as part of the job and automatically marked as running), and `PUT /jobs/<type>/<id>` to change its state.
