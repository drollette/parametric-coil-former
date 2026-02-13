# W7HAK Coil Former – WASM Plugin

A parametric coil-former generator for ham radio phasing lines, built in **Rust** and compiled to **WebAssembly**. Designed to run on **Cloudflare Pages** and be embedded as a plugin on [w7hak.com](https://w7hak.com).

## Features

- **Real-time 3D preview** via Three.js – rotate, zoom, and pan the model
- **Interactive sliders** for all parameters (wire length, wire diameter, PVC ID, pitch, rib clearance)
- **STL export** – download the model directly from the browser for 3D printing
- **Runs entirely client-side** – the WASM engine does all computation in the browser
- **Embeddable** – drop a single `<script>` tag into any page

## Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| Wire Length | 668 mm | 50–2000 mm | Target electrical phase length |
| Wire Diameter | 3.2 mm | 0.5–6.0 mm | Wire/coax core diameter |
| PVC Inner Diameter | 22.3 mm | 10–80 mm | Friction-fit target pipe ID |
| Pitch | 8.9 mm | 2–25 mm | Vertical spacing between wraps |
| Rib Clearance | 8.0 mm | 2–20 mm | Solid space at top/bottom |

## Quick Start

### Prerequisites

- [Rust](https://rustup.rs/) (stable, with `wasm32-unknown-unknown` target)
- [wasm-pack](https://rustwasm.github.io/wasm-pack/installer/)
- [Node.js](https://nodejs.org/) (for Wrangler CLI)

### Build

```bash
cd wasm-coil-former
make build
```

### Local Development

```bash
make dev
# Opens at http://localhost:8788
```

### Run Tests

```bash
make test
```

### Deploy to Cloudflare Pages

```bash
# Set your Cloudflare credentials first:
# export CLOUDFLARE_API_TOKEN=...
# export CLOUDFLARE_ACCOUNT_ID=...
make deploy
```

## Embedding on w7hak.com

Add the plugin to any page with:

```html
<div id="coil-former"></div>
<script src="https://coil-former.pages.dev/embed.js"
        data-target="#coil-former"
        data-height="700px"></script>
```

## Versioning

This project uses [Semantic Versioning](https://semver.org/). Releases are managed via GitHub Actions:

1. **CI** runs on every push/PR – checks, tests, and builds the WASM artifact
2. **Version Bump** – trigger manually via Actions to bump major/minor/patch
3. **Release** – pushing a `v*` tag builds, creates a GitHub Release with artifacts, and deploys to Cloudflare

To cut a new release:

```bash
# Via GitHub Actions UI: Actions → Version Bump → Run workflow
# Or manually:
# 1. Update version in Cargo.toml and package.json
# 2. Commit: git commit -am "chore: bump version to 0.2.0"
# 3. Tag:    git tag v0.2.0
# 4. Push:   git push origin main --tags
```

## Project Structure

```
wasm-coil-former/
├── Cargo.toml          # Rust package manifest
├── Makefile            # Build/dev/deploy shortcuts
├── package.json        # Node deps (wrangler)
├── wrangler.toml       # Cloudflare Workers config
├── rust-toolchain.toml # Pinned Rust toolchain
├── worker.js           # Cloudflare Worker entry (static serving)
├── src/
│   ├── lib.rs          # WASM entry point (wasm-bindgen exports)
│   └── geometry.rs     # Core geometry engine (math, mesh, STL)
└── static/
    ├── index.html      # Full-page app with Three.js viewport
    └── embed.js        # Lightweight iframe embedder for plugins
```

## License

MIT – see [LICENSE](../LICENSE)
