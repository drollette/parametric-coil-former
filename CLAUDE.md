# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Parametric generator for 3D-printable phasing coil formers (HAM radio antenna phasing lines). Python/CadQuery backend with a static HTML/JS frontend, containerized with Docker.

## Build & Run

```bash
# Docker (recommended)
docker-compose up --build
# Access: http://localhost:8000 (UI) | http://localhost:8000/docs (API docs)

# Local development
mamba env create -f environment.yml
mamba activate coil-former
pip install -r requirements.txt
uvicorn backend.main:app --reload --host 0.0.0.0 --port 8000

# Standalone CLI (generates STEP file to outputs/)
python phasing_coil.py
```

## Testing

No automated test suite exists. Test manually via the web UI or curl:
```bash
curl -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{"wire_len": 90.83, "wire_diam": 3.5, "pvc_id": 23.5}'
```

## Architecture

**Backend (Python/FastAPI):**
- `backend/main.py` — FastAPI app. `POST /generate` creates a UUID job, calls geometry engine, exports STL+STEP, schedules 1-hour cleanup. Serves frontend as static files from `wasm-coil-former/static/`.
- `backend/geometry.py` — CadQuery geometry engine. `build_coil_former()` creates the parametric 3D model (cylinder, V-groove helix, tunnels, optional friction ribs, chamfers). Returns a `CoilInfo` dataclass with computed values.
- `backend/schemas.py` — Pydantic models (`CoilRequest`, `CoilResponse`, `ComputedValues`) with validation ranges for all parameters.

**Frontend:**
- `wasm-coil-former/static/index.html` — Single-file SPA (vanilla HTML/CSS/JS, no framework). Dark theme UI with parameter sliders, 3D preview, and STL/STEP download buttons. Calls `POST /generate` on the backend.

**Standalone CLI:**
- `phasing_coil.py` — Self-contained script with hardcoded parameters at top of file. Generates STEP output directly. Predates the backend but still functional.

**Infrastructure:**
- `Dockerfile` — Based on `condaforge/mambaforge`, installs CadQuery via conda + FastAPI via pip.
- `docker-compose.yml` — Mounts `./outputs` and `./wasm-coil-former/static` (read-only, enables frontend hot reload).
- `environment.yml` — Conda env: Python 3.10 + CadQuery.
- `requirements.txt` — Pip: FastAPI, uvicorn, python-multipart.

## Key Design Details

- Coil diameter is auto-capped to `pvc_id - wire_diam` for PVC clearance safety.
- Turns calculated as: `sqrt(wire_len² / (circumference² + pitch²))`.
- CadQuery uses boolean `.cut()` operations to subtract grooves, tunnels, and bores from the base cylinder.
- The `wasm-coil-former/` directory name is historical — it previously held a Rust/WASM implementation. The frontend HTML is still served from there.

## CI/CD Note

The `.github/workflows/` CI and release workflows are outdated — they target the old Rust/WASM architecture (cargo, clippy, wasm-pack, Cloudflare Workers). They do not apply to the current Python/Docker stack.
