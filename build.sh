#!/usr/bin/bash
set -eu
# Build and push podman images in parallel.
user=jefftaylor42
images=${@:-dl align merge}
podman build -f docker/Dockerfile -t serratus-base .
for img in $images; do
    (
        podman build -f "docker/serratus-$img/Dockerfile" -t $user/serratus-$img .
        podman push $user/serratus-$img
    ) &
done
wait
