PROJECT := "rasusa"

# run clippy to check for linting issues
lint:
    cargo clippy --all-features --all-targets -- -D warnings

# run all tests
test:
    cargo test -v --all-targets --no-fail-fast

# get coverage with tarpaulin
coverage:
    cargo tarpaulin -t 300 -- --test-threads 1 --out xml