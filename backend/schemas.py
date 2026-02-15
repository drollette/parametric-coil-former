"""Pydantic request/response models for the Coil Former API."""

from pydantic import BaseModel, Field


class CoilRequest(BaseModel):
    """Request parameters for generating a coil former."""

    wire_len: float = Field(
        default=90.83,
        ge=5.0,
        le=2000.0,
        description="Target wire length in mm"
    )
    wire_diam: float = Field(
        default=3.5,
        ge=0.5,
        le=10.0,
        description="Wire diameter in mm"
    )
    pvc_id: float = Field(
        default=23.5,
        ge=10.0,
        le=100.0,
        description="PVC pipe inner diameter in mm"
    )
    coil_diameter: float = Field(
        default=15.0,
        ge=5.0,
        le=80.0,
        description="Desired coil diameter in mm (capped for PVC clearance)"
    )
    pitch: float = Field(
        default=10.0,
        ge=2.0,
        le=50.0,
        description="Vertical distance per turn in mm"
    )
    end_buffer: float = Field(
        default=10.0,
        ge=2.0,
        le=30.0,
        description="Space for ribs and transitions at top/bottom in mm"
    )
    enable_ribs: bool = Field(
        default=True,
        description="Enable friction ribs for PVC grip"
    )
    chamfer_size: float = Field(
        default=0.5,
        ge=0.0,
        le=3.0,
        description="Chamfer size on top/bottom edges in mm"
    )
    tunnel_tol: float = Field(
        default=0.2,
        ge=0.0,
        le=1.0,
        description="Extra diameter clearance for wire tunnels in mm"
    )


class ComputedValues(BaseModel):
    """Computed geometry values returned by the API."""

    turns: float = Field(description="Number of coil turns")
    total_height: float = Field(description="Total former height in mm")
    coil_diameter: float = Field(description="Actual coil diameter (may be capped)")
    winding_height: float = Field(description="Height of the winding section in mm")


class CoilResponse(BaseModel):
    """Response containing generated file paths and computed values."""

    uuid: str = Field(description="Unique job identifier")
    computed: ComputedValues
    stl_url: str = Field(description="URL to download STL file")
    step_url: str = Field(description="URL to download STEP file")
