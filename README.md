# Phasing Coil Former Generator

A parametric generator for 3D-printable phasing coil formers, designed for HAM radio antenna phasing lines. Features a web UI with real-time 3D preview and exports to both STL and STEP formats.

## What It Does

Generates a cylindrical coil former with a precision V-groove helix that guides wire placement. The former is sized to friction-fit inside standard PVC pipe, making it easy to create weatherproof phasing coils for antenna arrays.

## Features

- **Web Interface**: Real-time 3D preview with adjustable parameters
- **Dual Export**: Download STL for 3D printing or STEP for CAD editing
- **V-Groove Helix**: Precision groove sized to cradle wire securely
- **Straight Wire Tunnels**: Clean entry/exit paths from groove to center bore
- **Friction Ribs**: Optional ribs for secure fit inside PVC pipe
- **Safety Capping**: Coil diameter automatically limited for PVC clearance

## Quick Start (Docker)

```bash
# Build and run
docker build -t coil-former .
docker run -p 8000:8000 coil-former

# Or use docker-compose
docker-compose up --build
```

Then open http://localhost:8000

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `wire_len` | Target wire length in mm | 90.83 |
| `wire_diam` | Wire diameter in mm | 3.5 |
| `pvc_id` | PVC pipe inner diameter in mm | 23.5 |
| `coil_diameter` | Former diameter (auto-capped for clearance) | 15.0 |
| `pitch` | Vertical spacing between wraps in mm | 10.0 |
| `end_buffer` | Space for ribs/transitions at ends in mm | 10.0 |
| `chamfer_size` | Edge chamfer size in mm | 0.5 |
| `enable_ribs` | Add friction ribs for PVC grip | true |

## Local Python Usage

For command-line generation without the web UI:

```bash
# Set up environment
mamba env create -f environment.yml
mamba activate coil-former

# Edit parameters in phasing_coil.py, then run
python phasing_coil.py
```

Output files are saved to the `outputs/` directory.

## Project Structure

```
├── backend/
│   ├── main.py          # FastAPI application
│   ├── geometry.py      # CadQuery geometry engine
│   └── schemas.py       # Pydantic models
├── wasm-coil-former/
│   └── static/
│       └── index.html   # Web frontend
├── phasing_coil.py      # Standalone CLI script
├── Dockerfile
├── docker-compose.yml
├── environment.yml      # Conda dependencies
└── requirements.txt     # Python dependencies
```

## API

**POST /generate**

```json
{
  "wire_len": 90.83,
  "wire_diam": 3.5,
  "pvc_id": 23.5,
  "coil_diameter": 15.0,
  "pitch": 10.0,
  "end_buffer": 10.0,
  "enable_ribs": true,
  "chamfer_size": 0.5
}
```

Returns URLs to download generated STL and STEP files.

## License

MIT License - See [LICENSE](LICENSE) file

## Author

W7HAK - https://w7hak.com
