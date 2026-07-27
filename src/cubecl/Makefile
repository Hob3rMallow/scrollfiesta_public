SHELL := /bin/bash

.PHONY: fmt test rust-cpu cube-cpu hip bench-rust bench-cube-cpu bench-hip symbol clean

fmt:
	cargo fmt --all -- --check

# Fast correctness baseline; does not require CubeCL or ROCm.
test:
	cargo test --release --features rust-cpu

rust-cpu:
	cargo build --release --features rust-cpu

cube-cpu:
	cargo build --release --no-default-features --features rust-cpu,cubecl-cpu

hip:
	cargo build --release --no-default-features --features rust-cpu,cubecl-hip

bench-rust:
	cargo run --release --features rust-cpu --bin mls-bench -- rust-cpu

bench-cube-cpu:
	cargo run --release --no-default-features --features rust-cpu,cubecl-cpu --bin mls-bench -- cubecl-cpu

bench-hip:
	cargo run --release --no-default-features --features rust-cpu,cubecl-hip --bin mls-bench -- cubecl-hip

symbol: rust-cpu
	nm -D target/release/libherculaneum_mls_cubecl.so | grep ' T MLS_project_verts$$'

clean:
	cargo clean
