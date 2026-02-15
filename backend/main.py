"""
FastAPI application for the Coil Former generator.
Serves both the API and static frontend.
"""

import asyncio
import shutil
import uuid
from pathlib import Path

from fastapi import FastAPI, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .schemas import CoilRequest, CoilResponse, ComputedValues
from .geometry import build_coil_former, export_step, export_stl

# Paths
BASE_DIR = Path(__file__).parent.parent
OUTPUTS_DIR = BASE_DIR / "outputs"
STATIC_DIR = BASE_DIR / "wasm-coil-former" / "static"

# Ensure outputs directory exists
OUTPUTS_DIR.mkdir(exist_ok=True)

app = FastAPI(
    title="Coil Former API",
    description="Parametric coil former generator with CadQuery",
    version="2.0.0",
)

# CORS for development (allows localhost origins)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


async def cleanup_job(job_id: str, delay: int = 3600) -> None:
    """Delete job directory after specified delay (default 1 hour)."""
    await asyncio.sleep(delay)
    job_dir = OUTPUTS_DIR / job_id
    shutil.rmtree(job_dir, ignore_errors=True)


@app.post("/generate", response_model=CoilResponse)
async def generate_coil(
    request: CoilRequest,
    background_tasks: BackgroundTasks
) -> CoilResponse:
    """
    Generate a coil former with the specified parameters.
    Returns URLs to download STL and STEP files.
    """
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    job_dir = OUTPUTS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Build geometry
    result, info = build_coil_former(
        wire_len=request.wire_len,
        wire_diam=request.wire_diam,
        pvc_id=request.pvc_id,
        coil_diameter=request.coil_diameter,
        pitch=request.pitch,
        end_buffer=request.end_buffer,
        enable_ribs=request.enable_ribs,
        chamfer_size=request.chamfer_size,
        tunnel_tol=request.tunnel_tol,
    )

    # Export files
    stl_path = job_dir / "coil.stl"
    step_path = job_dir / "coil.step"

    export_stl(result, stl_path)
    export_step(result, step_path)

    # Schedule cleanup
    background_tasks.add_task(cleanup_job, job_id)

    return CoilResponse(
        uuid=job_id,
        computed=ComputedValues(
            turns=round(info.turns, 2),
            total_height=round(info.total_height, 2),
            coil_diameter=round(info.coil_diameter, 2),
            winding_height=round(info.winding_height, 2),
        ),
        stl_url=f"/outputs/{job_id}/coil.stl",
        step_url=f"/outputs/{job_id}/coil.step",
    )


@app.get("/outputs/{job_id}/{filename}")
async def download_file(job_id: str, filename: str) -> FileResponse:
    """Serve generated output files."""
    file_path = OUTPUTS_DIR / job_id / filename

    if not file_path.exists():
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="File not found")

    # Determine media type
    media_type = "application/octet-stream"
    if filename.endswith(".stl"):
        media_type = "model/stl"
    elif filename.endswith(".step"):
        media_type = "application/step"

    return FileResponse(
        path=file_path,
        media_type=media_type,
        filename=filename,
    )


@app.get("/health")
async def health_check() -> dict:
    """Health check endpoint."""
    return {"status": "healthy", "version": "2.0.0"}


# Serve static frontend files
# This must be mounted last to avoid catching API routes
app.mount("/", StaticFiles(directory=str(STATIC_DIR), html=True), name="static")
