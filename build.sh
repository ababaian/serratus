#!/bin/bash
set -eu
# Build and push container images in parallel.
images=${@:-dl align merge scheduler grafana prometheus}

DOCKERHUB_USER=${DOCKERHUB_USER:-}
DOCKER_BUILD=${DOCKER_BUILD:-sudo docker}

# Container Version
#VERSION=0 # Dev version
VERSION=0.1.2

# DOCKERHUB_USER='serratusbio'
# sudo docker login

$DOCKER_BUILD build -f docker/Dockerfile \
  -t serratus-base -t serratus-base:latest \
  -t serratus-base:$VERSION .

for img in $images; do
    (
        if [ "$img" = scheduler ]; then
            pushd scheduler
            $DOCKER_BUILD build  \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
            popd
        elif [ "$img" = grafana -o "$img" = prometheus ]; then
            pushd monitoring/"$img"
            $DOCKER_BUILD build  \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
            popd
        else
            $DOCKER_BUILD build -f "docker/serratus-$img/Dockerfile" \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
        fi

        if [ -z "$DOCKERHUB_USER" ]; then
          echo "No DOCKERHUB_USER set. Images are local only"
        else 
        # Push container images to repo
          $DOCKER_BUILD push $DOCKERHUB_USER/serratus-$img
          echo "Done pushing serratus-$img"
        fi
    ) &
done

wait
