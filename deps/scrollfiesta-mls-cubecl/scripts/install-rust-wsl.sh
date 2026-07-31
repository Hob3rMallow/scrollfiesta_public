#!/usr/bin/env bash
set -euo pipefail

if command -v cargo >/dev/null; then
    cargo --version
    exit 0
fi

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
# shellcheck disable=SC1090
source "$HOME/.cargo/env"
rustup toolchain install 1.88.0 --profile minimal
rustup default 1.88.0
cargo --version
