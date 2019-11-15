#!/usr/bin/env bash
set -eux

docker push "${IMAGE_NAME}:latest"
docker push "${IMAGE_NAME}:${TRAVIS_TAG}"
