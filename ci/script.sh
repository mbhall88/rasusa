# This script takes care of testing

set -ex

main() {
    cross fmt -- --check
    cross clippy --all-features --all-targets -- -D warnings
    cross build --target "$TARGET"
    cross build --target "$TARGET" --release

    if [ ! -z $DISABLE_TESTS ]; then
        return
    fi

    cross test --target "$TARGET" --all
    cross test --target "$TARGET" --release --all

    cross run --target "$TARGET" -- --help
    cross run --target "$TARGET" --release -- --help
}

# we don't run the "test phase" when doing deploys
if [ -z $TRAVIS_TAG ]; then
    main
fi
