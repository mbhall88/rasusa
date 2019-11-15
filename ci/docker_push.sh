#!/usr/bin/env bash
set -ex

docker push "${IMAGE_NAME}:latest"
docker push "${IMAGE_NAME}:${TRAVIS_TAG}"
