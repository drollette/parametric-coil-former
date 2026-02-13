"""
Phasing Coil Former Generator – CadQuery Port
Original OpenSCAD design by W7HAK

Generates a cylindrical coil former with a precision V-groove helix
for guiding wire placement.  Sized to friction-fit inside standard PVC pipe.
"""

import cadquery as cq
import math


def build_coil_former(
    wire_len=668.0,
    wire_diam=3.2,
    pvc_inner_diam=22.3,
    pitch=8.9,
    rib_clearance=8.0,
):
    """Build a parametric coil former and return the CadQuery solid.

    Parameters
    ----------
    wire_len : float
        Target wire length in mm (electrical phase length).
    wire_diam : float
        Wire diameter in mm (3.2 for RG58 core, 1.6 for bare wire).
    pvc_inner_diam : float
        Inside diameter of your PVC pipe in mm.
    pitch : float
        Vertical spacing between wraps in mm.
    rib_clearance : float
        Solid space at top & bottom for friction ribs in mm.

    Returns
    -------
    cq.Workplane
        The finished coil-former solid, ready for display or STL export.
    dict
        Computed dimensions (turns, coil height, total height, etc.).
    """

    # ── derived dimensions ───────────────────────────────────
    rib_diam      = pvc_inner_diam
    cylinder_diam = rib_diam - wire_diam
    cylinder_r    = cylinder_diam / 2.0
    center_bore_r = (wire_diam + 0.2) / 2.0
    r_wire        = wire_diam / 2.0

    circumference = math.pi * cylinder_diam
    calc_turns    = math.sqrt(wire_len ** 2 /
                              (circumference ** 2 + pitch ** 2))
    wh            = calc_turns * pitch            # winding height
    total_height  = wh + 2 * rib_clearance

    start_z   = rib_clearance
    end_z     = start_z + wh
    chamfer_h = (rib_diam - cylinder_diam) / 2.0
    v_depth   = r_wire * math.sqrt(2)             # 90° V cradles 50 %

    info = dict(
        turns=calc_turns,
        coil_height=wh,
        total_height=total_height,
    )

    # ── main body ────────────────────────────────────────────
    body = cq.Workplane("XY").circle(cylinder_r).extrude(total_height)

    # bottom friction rib (z = 2 mm)
    bot_z = 2.0
    bottom_rib = (
        cq.Workplane("XZ")
        .moveTo(cylinder_r, bot_z)
        .lineTo(rib_diam / 2, bot_z + chamfer_h)
        .lineTo(rib_diam / 2, bot_z + chamfer_h + 2)
        .lineTo(cylinder_r, bot_z + 2 * chamfer_h + 2)
        .close()
        .revolve(360, (0, 0), (0, 1))
    )

    # top friction rib (z = total_height − 8 mm)
    top_z = total_height - 8.0
    top_rib = (
        cq.Workplane("XZ")
        .moveTo(cylinder_r, top_z)
        .lineTo(rib_diam / 2, top_z + chamfer_h)
        .lineTo(rib_diam / 2, top_z + chamfer_h + 2)
        .lineTo(cylinder_r, top_z + 2 * chamfer_h + 2)
        .close()
        .revolve(360, (0, 0), (0, 1))
    )

    result = body.union(bottom_rib).union(top_rib)

    # centre bore
    bore = cq.Workplane("XY").circle(center_bore_r).extrude(total_height)
    result = result.cut(bore)

    # ── V-groove helix ───────────────────────────────────────
    helix = cq.Wire.makeHelix(pitch, wh, cylinder_r, lefthand=True)

    tan_y   = -2.0 * math.pi * cylinder_r
    tan_z   = pitch
    tan_len = math.hypot(tan_y, tan_z)
    tangent = cq.Vector(0, tan_y / tan_len, tan_z / tan_len)

    profile_plane = cq.Plane(
        origin=cq.Vector(cylinder_r, 0, 0),
        xDir=cq.Vector(1, 0, 0),
        normal=tangent,
    )

    oc = 1.0  # overshoot for clean boolean
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

    # ── entry / exit tunnels (purely radial for clean round holes) ──
    def _radial_tunnel(z, angle=0.0):
        overshoot = 2.0
        length = cylinder_r + center_bore_r + 2 * overshoot
        dir_vec = cq.Vector(math.cos(angle), math.sin(angle), 0)
        start_pt = cq.Vector(0, 0, z) - dir_vec * (center_bore_r + overshoot)
        solid = cq.Solid.makeCylinder(r_wire, length, pnt=start_pt, dir=dir_vec)
        return cq.Workplane("XY").newObject([solid])

    result = result.cut(_radial_tunnel(start_z, angle=0.0))

    end_angle = -2.0 * math.pi * calc_turns
    result = result.cut(_radial_tunnel(end_z, angle=end_angle))

    return result, info


# ── default build (for CQ-editor or ``from phasing_coil import result``) ──
result, info = build_coil_former()

print("=" * 40)
print(f"  Turns       : {info['turns']:.2f}")
print(f"  Coil height : {info['coil_height']:.2f} mm")
print(f"  Total height: {info['total_height']:.2f} mm")
print("=" * 40)
