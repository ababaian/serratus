#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${AWS_ACCESS_KEY_ID}" ]]; then
    echo 'Set AWS_ACCESS_KEY_ID'
    exit 1
fi

if [[ -z "${AWS_SECRET_ACCESS_KEY}" ]]; then
    echo 'Set AWS_SECRET_ACCESS_KEY'
    exit 1
fi

if [[ "${PWD##*/}" != "containers" ]]; then
    echo 'Must run from <repo_root>/containers'
    exit 1
fi

docker build -t serratus-align:test -f ./serratus-align/Dockerfile .

docker build -t serratus-align-runner -f ./serratus-align/test/Dockerfile.test .

docker run -ti \
    -e "AWS_ACCESS_KEY_ID" \
    -e "AWS_SECRET_ACCESS_KEY" \
    --name serratus-align-test \
    serratus-align-runner 

docker cp serratus-align-test:/home/serratus/UNIT_TEST.1337.bam ./serratus-align/test/TEST_RESULT.bam

docker rm serratus-align-test || true

diff ./serratus-align/test/EXPECTED.bam.golden \
    ./serratus-align/test/TEST_RESULT.bam
if [[ $? != 0 ]]; then
    echo 'Output differed from expected'
    echo 'See ./serratus-align/test/TEST_RESULT.bam for more info'
    exit 1
fi

echo 'Output matched expected!'

rm ./serratus-align/test/TEST_RESULT.bam
