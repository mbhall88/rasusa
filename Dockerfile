FROM rust:1.75 AS builder

COPY . /rasusa

WORKDIR /rasusa

ARG TARGET="x86_64-unknown-linux-musl"

RUN apt update \
    && apt install -y musl-tools \
    && rustup target add "$TARGET" \
    && cargo build --release --target "$TARGET" \
    && strip target/${TARGET}/release/rasusa


FROM bash:5.1

ARG TARGET="x86_64-unknown-linux-musl"
COPY --from=builder /rasusa/target/${TARGET}/release/rasusa /bin/

RUN rasusa --version

ENTRYPOINT [ "rasusa" ]

