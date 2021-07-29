# 1: Build the binary
FROM rust:1.53 as builder

# 1a: Prepare for static linking
RUN apt-get update && \
    apt-get dist-upgrade -y && \
    apt-get install -y musl-tools && \
    rustup target add x86_64-unknown-linux-musl

# 1b: Download and compile Rust dependencies (and store as a separate Docker layer)
COPY . /rasusa
WORKDIR /rasusa
RUN cargo install --target x86_64-unknown-linux-musl --path .

# 2: Copy the binary to an empty Docker image
FROM bash:5.0

COPY --from=builder /usr/local/cargo/bin/rasusa /bin/rasusa

RUN rasusa --help

CMD ["/bin/rasusa"]
