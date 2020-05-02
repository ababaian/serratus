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
