#!/usr/bin/env bash
# This script takes care of testing

set -eux

main() {
    cross fmt -- --check
    cross clippy --all-features --all-targets -- -D warnings
    cross build
    cross build --release

    if [ -n "$DISABLE_TESTS" ]; then
        return
    fi

    cross test --all
    cross test --release --all

    cross run -- --help
    cross run --release -- --help
}

# we don't run the "test phase" when doing deploys
if [ -z "$TRAVIS_TAG" ]; then
    main
fi
