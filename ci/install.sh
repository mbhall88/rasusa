#!/usr/bin/env sh
set -ex

main() {
    if [ "$TRAVIS_OS_NAME" = osx ]; then
        alias sort='gsort'  # for `sort --sort-version`, from brew's coreutils.
    fi

    if [ "$TRAVIS_RUST_VERSION" = nightly ]; then
        rustup toolchain remove stable && rustup toolchain add stable
        rustup component add clippy --toolchain=nightly || cargo install --git https://github.com/rust-lang/rust-clippy/ --force clippy
        rustup component add rustfmt --toolchain nightly
    else
        rustup component add clippy rustfmt
    fi
    # This fetches latest stable release of cross
    tag=$(git ls-remote --tags --refs --exit-code https://github.com/japaric/cross \
                       | cut -d/ -f3 \
                       | grep -E '^v[0.1.0-9.]+$' \
                       | sort --version-sort \
                       | tail -n1)
    curl -LSfs https://japaric.github.io/trust/install.sh | \
        sh -s -- \
           --force \
           --git japaric/cross \
           --tag "$tag" \
           --target "$TARGET"
}

main
