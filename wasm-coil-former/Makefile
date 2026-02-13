.PHONY: build test clean dev deploy fmt lint

# Build the WASM package and assemble dist/
build:
	wasm-pack build --target web --out-dir pkg --release
	mkdir -p dist
	cp static/* dist/
	cp pkg/wasm_coil_former_bg.wasm dist/
	cp pkg/wasm_coil_former.js dist/

# Run Rust unit tests
test:
	cargo test --verbose

# Lint & format checks
lint:
	cargo clippy --target wasm32-unknown-unknown -- -D warnings
	cargo fmt -- --check

# Auto-format
fmt:
	cargo fmt

# Local dev server via wrangler
dev: build
	npx wrangler pages dev dist

# Deploy to Cloudflare Pages
deploy: build
	npx wrangler pages deploy dist --project-name=coil-former

# Clean build artifacts
clean:
	cargo clean
	rm -rf pkg dist node_modules
