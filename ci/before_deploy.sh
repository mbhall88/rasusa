#!/usr/bin/env sh
# This script takes care of building the crate, packaging it for release, and building
# a docker image from the specified release binary

set -ex

build_docker() {
    DOCKER_BIN="target/docker"

    mkdir -p "$DOCKER_BIN"
    cp "target/${TARGET}/release/${PROJECT_NAME}" "${DOCKER_BIN}/${PROJECT_NAME}"

    docker build --tag "$IMAGE_NAME" .

    docker images

    echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin

    docker tag "$IMAGE_NAME" "${IMAGE_NAME}:latest"
    docker tag "$IMAGE_NAME" "${IMAGE_NAME}:${TRAVIS_TAG}"
}

main() {
    src=$(pwd)
    stage=

    case "$TRAVIS_OS_NAME" in
        linux)
            stage=$(mktemp -d)
            ;;
        osx)
            stage=$(mktemp -d -t tmp)
            ;;
    esac

    test -f Cargo.lock || cargo generate-lockfile

    cross build --target "$TARGET" --release

    if [ "$TARGET" = x86_64-unknown-linux-musl ]; then
        build_docker
    fi

    cp "target/${TARGET}/release/${PROJECT_NAME}" "${stage}/"

    cd "$stage" || exit 1
    tar czf "${src}/${PROJECT_NAME}-${TRAVIS_TAG}-${TARGET}.tar.gz" ./*
    cd "$src" || exit 1

    rm -rf "$stage"
}

main
