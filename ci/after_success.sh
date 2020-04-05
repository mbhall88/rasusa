#!/usr/bin/env bash
# This script takes care of uploading code coverage to codecov
set -ex

main() {
    kcov_version="38"
    url="https://github.com/SimonKagstrom/kcov/archive/${kcov_version}.tar.gz"
    wget "$url" -O - | tar xzf -
    cd kcov-${kcov_version} || exit 1
    mkdir build
    cd build || exit 1
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make || exit 1
    make install DESTDIR=../../kcov-build || exit 1
    cd ../../ || exit 1
    rm -rf kcov-${kcov_version}
    cargo install --force cargo-kcov
    cargo kcov -v --kcov ./kcov-build/usr/local/bin/kcov -- --verify --exclude-pattern=/.cargo,/usr/lib,main.rs
    bash <(curl -s https://codecov.io/bash) &&
    echo "Uploaded code coverage"
}

# we don't run tarpaulin on osx and only run on stable
if [[ "$TRAVIS_OS_NAME" == linux && "$TRAVIS_RUST_VERSION" == stable ]]; then
    main
fi
