#!/bin/bash
set -eu
# Build and push container images in parallel.
images=${@:-dl align merge scheduler grafana prometheus}
DOCKERHUB_USER=${DOCKERHUB_USER:-}

# Use `podman` or `docker` as builder
#builder='podman'
builder='sudo docker'

# Container Version
#VERSION=0 # Dev version
VERSION=0.1.2

# DOCKERHUB_USER='serratusbio'
# sudo docker login

if [ -z "$DOCKERHUB_USER" ]; then
    echo "Environment variable DOCKERHUB_USER not set."
    exit 1
fi

$builder build -f docker/Dockerfile \
  -t serratus-base -t serratus-base:latest \
  -t serratus-base:$VERSION .

for img in $images; do
    (
        if [ "$img" = scheduler ]; then
            pushd scheduler
            $builder build  \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
            popd
        elif [ "$img" = grafana -o "$img" = prometheus ]; then
            pushd monitoring/"$img"
            $builder build  \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
            popd
        else
            $builder build -f "docker/serratus-$img/Dockerfile" \
              -t serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img \
              -t $DOCKERHUB_USER/serratus-$img:$VERSION \
              -t $DOCKERHUB_USER/serratus-$img:latest .
        fi

        $builder push $DOCKERHUB_USER/serratus-$img
        echo "Done pushing serratus-$img"
    ) &
done

wait
