#!/usr/bin/env sh
# This script takes care of building your crate and packaging it for release

set -ex

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

    cp "target/${TARGET}/release/${PROJECT_NAME}" "${stage}/"

    cd "$stage" || exit 1
    tar czf "${src}/${PROJECT_NAME}-${TRAVIS_TAG}-${TARGET}.tar.gz" ./*
    cd "$src" || exit 1

    rm -rf "$stage"
}

main
