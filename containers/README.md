# Serratus Dockerfile README

[ See: CONTRIBUTING.md for detailed documention ](https://github.com/ababaian/serratus/blob/master/CONTRIBUTING.md#production-containers-and-code)

## Quick References
```
# Development Notes:
# serratus-base
# |    |--serratus-dl
# |    |--serratus-align
# |    |--serratus-merge
# serratus-scheduler
# serratus-grafana
# serratus-prometheus
```

## Local pipeline
Serratus provides a simplified image for local development. This is a good place to start.

First, build the `serratus-batch` _image_. From this directory:
`docker build -t serratus-batch -f ./serratus-batch/Dockerfile .`
* `-t serratus-batch` - Tag the image so we can reference it later
* `-f ./serratus-batch/Dockerfile` - Use the indicated Dockerfile to build the image
* `.` - Root the build at the current directory

Now that you've got the _image_ built, you can run a new _container_ from that image. `serratus-batch` is a tools container, which means it is expected to run, and then exit. It writes output data to AWS s3.

For s3 access, you'll need to set up an IAM keypair. Set the following environment variables - they'll be passed to the running _container_.
```
export AWS_ACCESS_KEY_ID="AKIA...."
export AWS_SECRET_ACCESS_KEY="secret_key_here"
```

Then, run the _container_
`docker run --rm -ti -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY serratus-batch run <SRA> <Index>`
* `--rm` - Remove the container after it exits
* `-ti` - Hook up a TTY and Interactive mode. tl;dr - Get output on your screen, kill it with Ctrl+C
* `-e ...` - Forward your AWS credentials to the container
* `serratus-batch` - This is the image to use for the new container. Reference by the tag we set when we built the image
* `run <SRA> <Index>` - This is the command that will be executed in the container

Tips:
* To poke around the container, substitute `/bin/bash` for `run ...`. That will run the container and drop you into a bash shell, where you can run the pipeline piece-by-piece yourself.
* To save a pre-fetched accession for experimentation:
  1. Run with `/bin/bash`
  2. `prefetch <SRA>`
  3. From another shell, `docker ps` to find your container name. Then `docker commit <name> docker-batch:<pick_a_tag_name>`
  4. Substitute `docker-batch:tag` as the image for `docker run`, and `/bin/bash` as the command
  5. Manually run the pipeline. Copy commands from `serratus-batch/run` to start out.
* To build a tagged image with a pre-downloaded index, follow the above steps. Instead of `prefetch <SRA>`, run `aws s3 cp ...` (copy from `run.sh`), and commit the result

## Build all containers for serratus
Instructions for EC2/Amazon Linux
```
# From base amazon linux 2
sudo yum install -y docker
sudo yum install -y git
sudo service docker start
```

```
git clone https://github.com/ababaian/serratus.git; cd serratus

# If you want to upload containers to your repository, include this.
export DOCKERHUB_USER='serratusbio' # optional
sudo docker login # optional

# Build all containers and upload them docker hub repo (if available)
./build_containers.sh   # run this in the folder 'serratus/containers'

```

## Run interactive serratus-dl
```
sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:latest
```
