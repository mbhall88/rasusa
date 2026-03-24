# External Integrations

**Analysis Date:** 2025-01-24

## APIs & External Services

**CI/CD:**
- GitHub Actions - Automates testing, building, and publishing.
  - Auth: `GITHUB_TOKEN` (automatically provided).

## Data Storage

**Databases:**
- None detected.

**File Storage:**
- Local filesystem only - Handles high-throughput reading and writing of sequence (FASTQ/FASTA) and alignment (SAM/BAM/CRAM) data.

**Caching:**
- GitHub Actions cache - Used for caching Cargo dependencies and build artifacts to speed up CI.

## Authentication & Identity

**Auth Provider:**
- Custom (Token-based)
  - `GITHUB_TOKEN` - Used for GitHub Packages and Releases.
  - `CARGO_REGISTRY_TOKEN` - Used for publishing to `crates.io`.
  - `CODECOV_TOKEN` - Used for uploading coverage reports to `codecov.io`.

## Monitoring & Observability

**Error Tracking:**
- None (standard stderr output).

**Logs:**
- Standard log output via `env_logger`.

## CI/CD & Deployment

**Hosting:**
- `ghcr.io` (GitHub Container Registry) - Hosting Docker images.
- `crates.io` - Hosting Rust packages.
- GitHub Releases - Hosting compiled binaries for various architectures.

**CI Pipeline:**
- GitHub Actions (`.github/workflows/rust-ci.yaml`, `.github/workflows/docker.yaml`, `.github/workflows/release.yaml`)
- `release-please` - Automated release management and versioning.
- `codecov` - Integration for tracking code coverage.

## Environment Configuration

**Required env vars:**
- `GITHUB_TOKEN` (CI only)
- `CARGO_REGISTRY_TOKEN` (CI only)
- `CODECOV_TOKEN` (CI only)
- `RUST_LOG` (Runtime) - Optional, for controlling log level.

**Secrets location:**
- GitHub Repository Secrets.

## Webhooks & Callbacks

**Incoming:**
- GitHub repository events (push, pull_request, release).

**Outgoing:**
- None.

---

*Integration audit: 2025-01-24*
