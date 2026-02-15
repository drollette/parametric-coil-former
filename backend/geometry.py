"""
CadQuery geometry engine for generating coil formers.
Refactored from phasing_coil.py for use as a library.
"""

import cadquery as cq
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CoilInfo:
    """Computed geometry values for a coil former."""
    turns: float
    total_height: float
    coil_diameter: float
    winding_height: float


def build_coil_former(
    wire_len: float = 90.83,
    wire_diam: float = 3.5,
    pvc_id: float = 23.5,
    coil_diameter: float = 15.0,
    pitch: float = 10.0,
    end_buffer: float = 10.0,
    enable_ribs: bool = True,
    chamfer_size: float = 0.5,
    tunnel_tol: float = 0.2,
) -> tuple[cq.Workplane, CoilInfo]:
    """
    Build a parametric coil former with V-groove helix and wire tunnels.

    Args:
        wire_len: Target wire length in mm
        wire_diam: Wire diameter in mm
        pvc_id: PVC pipe inner diameter in mm
        coil_diameter: Desired coil diameter (capped for clearance)
        pitch: Vertical distance per turn in mm
        end_buffer: Space for ribs/transitions at ends in mm
        enable_ribs: Whether to add friction ribs
        chamfer_size: Chamfer on top/bottom edges in mm
        tunnel_tol: Extra clearance for wire tunnels in mm

    Returns:
        Tuple of (CadQuery Workplane result, CoilInfo with computed values)
    """
    # Safety validation - cap coil diameter for PVC clearance
    max_safe_diam = pvc_id - wire_diam
    actual_coil_diam = min(coil_diameter, max_safe_diam)

    # Derived dimensions
    cylinder_r = actual_coil_diam / 2.0
    rib_diam = pvc_id
    center_bore_r = (wire_diam + 0.2) / 2.0
    r_wire = wire_diam / 2.0
    r_tunnel = r_wire + (tunnel_tol / 2.0)

    # Calculate turns from wire length
    circumference = math.pi * actual_coil_diam
    base_turns = math.sqrt(wire_len ** 2 / (circumference ** 2 + pitch ** 2))
    extension_turns = wire_diam / circumference
    total_turns = base_turns + extension_turns

    # Heights
    winding_height = total_turns * pitch
    total_height = winding_height + (2 * end_buffer)

    start_z = end_buffer
    end_z = start_z + winding_height
    v_depth = r_wire * math.sqrt(2)

    # Build main body cylinder
    result = cq.Workplane("XY").circle(cylinder_r).extrude(total_height)

    # Chamfer top/bottom edges
    if chamfer_size > 0:
        result = result.edges(">Z or <Z").chamfer(chamfer_size)

    # Friction ribs (conditional)
    if enable_ribs:
        def _add_rib(z_pos: float) -> cq.Workplane:
            rib_height_diff = (rib_diam - actual_coil_diam) / 2.0
            return (
                cq.Workplane("XZ")
                .moveTo(cylinder_r, z_pos)
                .lineTo(rib_diam / 2, z_pos + rib_height_diff)
                .lineTo(rib_diam / 2, z_pos + rib_height_diff + 2)
                .lineTo(cylinder_r, z_pos + 2 * rib_height_diff + 2)
                .close()
                .revolve(360, (0, 0), (0, 1))
            )

        result = result.union(_add_rib(2.0))
        result = result.union(
            _add_rib(total_height - 6.0 - ((rib_diam - actual_coil_diam) / 2.0))
        )

    # V-groove helix
    helix = cq.Wire.makeHelix(pitch, winding_height, cylinder_r, lefthand=True)

    tan_y, tan_z = -2.0 * math.pi * cylinder_r, pitch
    tan_len = math.hypot(tan_y, tan_z)
    profile_plane = cq.Plane(
        origin=cq.Vector(cylinder_r, 0, 0),
        xDir=cq.Vector(1, 0, 0),
        normal=cq.Vector(0, tan_y / tan_len, tan_z / tan_len),
    )

    oc = 1.0
    groove = (
        cq.Workplane(profile_plane)
        .moveTo(-v_depth, 0)
        .lineTo(oc, v_depth + oc)
        .lineTo(oc, -(v_depth + oc))
        .close()
        .sweep(helix, isFrenet=True)
        .translate(cq.Vector(0, 0, start_z))
    )
    result = result.cut(groove)

    # Straight entry/exit tunnels
    def _straight_tunnel(z_level: float, angle: float, is_top: bool = False) -> cq.Workplane:
        """Create a straight tunnel from helix entry point to center bore."""
        z_internal_buffer = 2.0
        z_exit = (total_height - z_internal_buffer) if is_top else z_internal_buffer

        p1 = cq.Vector(
            cylinder_r * math.cos(angle),
            cylinder_r * math.sin(angle),
            z_level
        )
        p2 = cq.Vector(0, 0, z_exit)

        path = cq.Wire.assembleEdges([cq.Edge.makeLine(p1, p2)])
        tunnel_plane = cq.Plane(origin=p1, normal=p1.sub(p2))

        return (
            cq.Workplane(tunnel_plane)
            .circle(r_tunnel)
            .sweep(cq.Workplane().newObject([path]))
        )

    angle_offset = (extension_turns / 2.0) * 2.0 * math.pi

    # Bottom entry tunnel
    result = result.cut(
        _straight_tunnel(start_z + (wire_diam / 2), -angle_offset, False)
    )
    # Top exit tunnel
    result = result.cut(
        _straight_tunnel(
            end_z - (wire_diam / 2),
            (-2.0 * math.pi * total_turns) + angle_offset,
            True
        )
    )

    # Center bore
    bore = cq.Workplane("XY").circle(center_bore_r).extrude(total_height)
    result = result.cut(bore)

    info = CoilInfo(
        turns=total_turns,
        total_height=total_height,
        coil_diameter=actual_coil_diam,
        winding_height=winding_height,
    )

    return result, info


def export_step(result: cq.Workplane, output_path: Path) -> None:
    """Export CadQuery result to STEP format."""
    cq.exporters.export(result, str(output_path))


def export_stl(result: cq.Workplane, output_path: Path) -> None:
    """Export CadQuery result to STL format."""
    cq.exporters.export(result, str(output_path), exportType="STL")
