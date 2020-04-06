#!/usr/bin/bash
set -eu
# Build and push podman images in parallel.
images=${@:-dl align merge scheduler grafana prometheus}
DOCKERHUB_USER=${DOCKERHUB_USER:-}

if [ -z "$DOCKERHUB_USER" ]; then
    echo "Environment variable DOCKERHUB_USER not set."
    exit 1
fi

podman build -f docker/Dockerfile -t serratus-base .
for img in $images; do
    (
        if [ "$img" = scheduler ]; then
            pushd scheduler
            podman build -t $DOCKERHUB_USER/serratus-$img .
            popd
        elif [ "$img" = grafana -o "$img" = prometheus ]; then
            pushd monitoring/"$img"
            podman build -t $DOCKERHUB_USER/serratus-$img .
            popd
        else
            podman build -f "docker/serratus-$img/Dockerfile" -t $DOCKERHUB_USER/serratus-$img .
        fi

        podman push -q $DOCKERHUB_USER/serratus-$img
        echo "Done pushing serratus-$img"
    ) &
done

wait
