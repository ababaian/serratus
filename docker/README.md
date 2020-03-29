# Serratus Dockerfile README

```
# Development Notes:
# serratus-base
# |    |--serratus-dl
# |    |--serratus-align
# |    |     |--serratus-align-hgr1
# |    |     |--...
# |    |--serratus-merge
```

## Build all containers for serratus

```
git clone https://github.com/ababaian/serratus.git; cd serratus
sudo docker build -t serratus-base:0 -t serratus-base:latest -f docker/Dockerfile .
sudo docker build -t serratus-dl:0 -t serratus-dl:latest -f docker/serratus-dl/Dockerfile .
sudo docker build -t serratus-align:0 -t serratus-align:latest -f docker/serratus-align/Dockerfile .
sudo docker build -t serratus-merge:0 -t serratus-merge:latest -f docker/serratus-merge/Dockerfile .
```

## Uploading serratus container images to AWS ECR
Paste resulting command in terminal to authenticate

```
## NOTE: Depricated. Jeff builds to dockerhub directly
##       these instructions are for 
#aws ecr get-login --region us-east-1 --no-include-email

## Tag images with ECR information
# sudo docker tag serratus-base:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker tag serratus-dl:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker tag serratus-align:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker tag serratus-merge:latest 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge

## Push images to ECR (rather host on dockerhub)
## Make sure to initialize each repository on your AWS ECR account
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-base
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-dl
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-align
# sudo docker push 797308887321.dkr.ecr.us-east-1.amazonaws.com/serratus-merge
```

## Run interactive serratus-dl
```
sudo docker run --rm --entrypoint /bin/bash -it serratus-dl:0
```

## Testing scheduler
```
cd serratus/scheduler
sudo docker build -t scheduler:0 .
docker run -d --rm -p8000:8000 --name sch scheduler:0
curl -T /path/to/SraRunInfo.csv localhost:8000/add_sra_run_info
```
