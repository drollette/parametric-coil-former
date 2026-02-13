//! Core coil-former geometry engine.
//!
//! Translates the parametric OpenSCAD / CadQuery model into pure Rust.
//! All units are millimetres.

use std::f64::consts::PI;

/// User-facing parameters that fully define a coil former.
#[derive(Debug, Clone, Copy)]
pub struct CoilParams {
    /// Target wire length in mm (electrical phase length).
    pub wire_len: f64,
    /// Wire diameter in mm (3.2 for RG58 core, 1.6 for bare wire).
    pub wire_diam: f64,
    /// Inside diameter of the PVC pipe in mm.
    pub pvc_inner_diam: f64,
    /// Vertical spacing between wraps in mm.
    pub pitch: f64,
    /// Solid space at top & bottom for friction ribs in mm.
    pub rib_clearance: f64,
}

impl Default for CoilParams {
    fn default() -> Self {
        Self {
            wire_len: 668.0,
            wire_diam: 3.2,
            pvc_inner_diam: 22.3,
            pitch: 8.9,
            rib_clearance: 8.0,
        }
    }
}

/// Computed dimensions derived from [`CoilParams`].
#[derive(Debug, Clone, Copy)]
pub struct CoilDerived {
    pub rib_diam: f64,
    pub cylinder_diam: f64,
    pub cylinder_r: f64,
    pub center_bore_r: f64,
    pub r_wire: f64,
    pub circumference: f64,
    pub calc_turns: f64,
    pub winding_height: f64,
    pub total_height: f64,
    pub start_z: f64,
    pub end_z: f64,
    pub chamfer_h: f64,
    pub v_depth: f64,
}

impl CoilDerived {
    pub fn from_params(p: &CoilParams) -> Self {
        let rib_diam = p.pvc_inner_diam;
        let cylinder_diam = rib_diam - p.wire_diam;
        let cylinder_r = cylinder_diam / 2.0;
        let center_bore_r = (p.wire_diam + 0.2) / 2.0;
        let r_wire = p.wire_diam / 2.0;

        let circumference = PI * cylinder_diam;
        let calc_turns =
            (p.wire_len * p.wire_len / (circumference * circumference + p.pitch * p.pitch)).sqrt();
        let winding_height = calc_turns * p.pitch;
        let total_height = winding_height + 2.0 * p.rib_clearance;

        let start_z = p.rib_clearance;
        let end_z = start_z + winding_height;
        let chamfer_h = (rib_diam - cylinder_diam) / 2.0;
        let v_depth = r_wire * std::f64::consts::SQRT_2;

        Self {
            rib_diam,
            cylinder_diam,
            cylinder_r,
            center_bore_r,
            r_wire,
            circumference,
            calc_turns,
            winding_height,
            total_height,
            start_z,
            end_z,
            chamfer_h,
            v_depth,
        }
    }
}

// ─── Mesh primitives ───────────────────────────────────────────────

/// A simple triangle-soup mesh (position-only).
#[derive(Debug, Clone)]
pub struct TriMesh {
    /// Flat array of vertex positions: [x0,y0,z0, x1,y1,z1, …]
    pub positions: Vec<f32>,
    /// Triangle indices into the positions array (every 3 values = 1 tri).
    pub indices: Vec<u32>,
}

impl TriMesh {
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            indices: Vec::new(),
        }
    }

    /// Merge another mesh into this one.
    pub fn append(&mut self, other: &TriMesh) {
        let base = (self.positions.len() / 3) as u32;
        self.positions.extend_from_slice(&other.positions);
        self.indices
            .extend(other.indices.iter().map(|i| i + base));
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len() / 3
    }
}

// ─── Cylinder generation ───────────────────────────────────────────

/// Closed cylinder along the Z axis centred at origin.
pub fn make_cylinder(radius: f64, height: f64, segments: u32) -> TriMesh {
    let mut mesh = TriMesh::new();
    let segs = segments as usize;

    // --- vertices ---
    // bottom ring: 0..segs
    // top ring:    segs..2*segs
    // bottom cap centre: 2*segs
    // top cap centre:    2*segs+1
    for ring in 0..2 {
        let z = if ring == 0 { 0.0 } else { height };
        for i in 0..segs {
            let angle = 2.0 * PI * (i as f64) / (segs as f64);
            mesh.positions.push((radius * angle.cos()) as f32);
            mesh.positions.push((radius * angle.sin()) as f32);
            mesh.positions.push(z as f32);
        }
    }
    // cap centres
    mesh.positions.extend_from_slice(&[0.0, 0.0, 0.0]);
    mesh.positions.extend_from_slice(&[0.0, 0.0, height as f32]);

    let bot_center = (2 * segs) as u32;
    let top_center = (2 * segs + 1) as u32;

    for i in 0..segs {
        let next = (i + 1) % segs;
        let b0 = i as u32;
        let b1 = next as u32;
        let t0 = (segs + i) as u32;
        let t1 = (segs + next) as u32;

        // side quad (2 tris)
        mesh.indices.extend_from_slice(&[b0, t0, b1]);
        mesh.indices.extend_from_slice(&[b1, t0, t1]);

        // bottom cap (winding: CW from below → CCW from outside)
        mesh.indices.extend_from_slice(&[bot_center, b1, b0]);

        // top cap
        mesh.indices.extend_from_slice(&[top_center, t0, t1]);
    }

    mesh
}

/// Conical frustum (truncated cone) along Z axis.
pub fn make_frustum(r_bottom: f64, r_top: f64, height: f64, segments: u32) -> TriMesh {
    let mut mesh = TriMesh::new();
    let segs = segments as usize;

    for ring in 0..2 {
        let (r, z) = if ring == 0 {
            (r_bottom, 0.0)
        } else {
            (r_top, height)
        };
        for i in 0..segs {
            let angle = 2.0 * PI * (i as f64) / (segs as f64);
            mesh.positions.push((r * angle.cos()) as f32);
            mesh.positions.push((r * angle.sin()) as f32);
            mesh.positions.push(z as f32);
        }
    }
    mesh.positions.extend_from_slice(&[0.0, 0.0, 0.0]);
    mesh.positions
        .extend_from_slice(&[0.0, 0.0, height as f32]);

    let bot_center = (2 * segs) as u32;
    let top_center = (2 * segs + 1) as u32;

    for i in 0..segs {
        let next = (i + 1) % segs;
        let b0 = i as u32;
        let b1 = next as u32;
        let t0 = (segs + i) as u32;
        let t1 = (segs + next) as u32;

        mesh.indices.extend_from_slice(&[b0, t0, b1]);
        mesh.indices.extend_from_slice(&[b1, t0, t1]);
        mesh.indices.extend_from_slice(&[bot_center, b1, b0]);
        mesh.indices.extend_from_slice(&[top_center, t0, t1]);
    }

    mesh
}

// ─── Transform helpers ─────────────────────────────────────────────

pub fn translate(mesh: &mut TriMesh, dx: f32, dy: f32, dz: f32) {
    for i in (0..mesh.positions.len()).step_by(3) {
        mesh.positions[i] += dx;
        mesh.positions[i + 1] += dy;
        mesh.positions[i + 2] += dz;
    }
}

pub fn rotate_z(mesh: &mut TriMesh, angle_rad: f32) {
    let (s, c) = (angle_rad.sin(), angle_rad.cos());
    for i in (0..mesh.positions.len()).step_by(3) {
        let x = mesh.positions[i];
        let y = mesh.positions[i + 1];
        mesh.positions[i] = x * c - y * s;
        mesh.positions[i + 1] = x * s + y * c;
    }
}

// ─── Full coil former assembly ─────────────────────────────────────

const SEGMENTS: u32 = 64;
const HELIX_STEPS_PER_TURN: u32 = 60;

/// Build the complete coil former mesh from parameters.
pub fn build_coil_former(params: &CoilParams) -> (TriMesh, CoilDerived) {
    let d = CoilDerived::from_params(params);

    let mut mesh = TriMesh::new();

    // 1. Main body cylinder
    let body = make_cylinder(d.cylinder_r, d.total_height, SEGMENTS);
    mesh.append(&body);

    // 2. Friction ribs
    let bottom_rib = build_friction_rib(d.cylinder_r, d.rib_diam / 2.0, d.chamfer_h, 2.0);
    let top_rib = build_friction_rib(
        d.cylinder_r,
        d.rib_diam / 2.0,
        d.chamfer_h,
        d.total_height - 8.0,
    );
    mesh.append(&bottom_rib);
    mesh.append(&top_rib);

    // 3. Center bore (rendered as a thin inner cylinder for visual representation)
    let bore = make_cylinder(d.center_bore_r, d.total_height, SEGMENTS);
    mesh.append(&bore);

    // 4. V-groove helix (rendered as tube following helical path)
    let groove = build_helix_groove(&d);
    mesh.append(&groove);

    // 5. Entry/exit tunnels
    let entry = build_tunnel(d.cylinder_r, d.r_wire, d.start_z, 0.0);
    let exit_angle = -2.0 * PI * d.calc_turns;
    let exit = build_tunnel(d.cylinder_r, d.r_wire, d.end_z, exit_angle);
    mesh.append(&entry);
    mesh.append(&exit);

    (mesh, d)
}

/// Build a single friction rib (chamfer-up, flat, chamfer-down).
fn build_friction_rib(cyl_r: f64, rib_r: f64, chamfer_h: f64, base_z: f64) -> TriMesh {
    let mut mesh = TriMesh::new();

    // chamfer up
    let mut c1 = make_frustum(cyl_r, rib_r, chamfer_h, SEGMENTS);
    translate(&mut c1, 0.0, 0.0, base_z as f32);
    mesh.append(&c1);

    // flat band
    let mut flat = make_cylinder(rib_r, 2.0, SEGMENTS);
    translate(&mut flat, 0.0, 0.0, (base_z + chamfer_h) as f32);
    mesh.append(&flat);

    // chamfer down
    let mut c2 = make_frustum(rib_r, cyl_r, chamfer_h, SEGMENTS);
    translate(&mut c2, 0.0, 0.0, (base_z + chamfer_h + 2.0) as f32);
    mesh.append(&c2);

    mesh
}

/// Build the V-groove helix as a sequence of small oriented wedges along the
/// helical path.  For the WASM preview we approximate with small cones.
fn build_helix_groove(d: &CoilDerived) -> TriMesh {
    let total_steps = (d.calc_turns * HELIX_STEPS_PER_TURN as f64).round() as u32;
    let total_steps = total_steps.max(1);
    let dz = d.winding_height / total_steps as f64;
    let da = -2.0 * PI * d.calc_turns / total_steps as f64;

    let mut mesh = TriMesh::new();

    // Each segment: a small cone positioned on the helix path cut into the surface.
    let seg_len = {
        let arc = d.cylinder_r * da.abs();
        (arc * arc + dz * dz).sqrt()
    };

    for i in 0..total_steps {
        let z0 = d.start_z + (i as f64) * dz;
        let a0 = (i as f64) * da;
        let z1 = d.start_z + ((i + 1) as f64) * dz;
        let a1 = ((i + 1) as f64) * da;

        // Mid-point on helix surface
        let mid_z = ((z0 + z1) / 2.0) as f32;
        let mid_a = ((a0 + a1) / 2.0) as f32;

        // Small groove marker: a thin frustum pointing inward
        let mut seg = make_frustum(0.01, d.v_depth + 0.5, seg_len, 8);

        // Orient along the helix tangent (approximate: rotate about Y then Z)
        let pitch_angle = (dz / (d.cylinder_r * da.abs())).atan() as f32;
        rotate_y(&mut seg, -pitch_angle);

        // Position on cylinder surface
        translate(
            &mut seg,
            (d.cylinder_r - d.v_depth) as f32,
            0.0,
            0.0,
        );
        rotate_z(&mut seg, mid_a);
        translate(&mut seg, 0.0, 0.0, mid_z);

        mesh.append(&seg);
    }

    mesh
}

/// Build a radial tunnel cylinder for wire entry/exit.
fn build_tunnel(cyl_r: f64, r_wire: f64, z: f64, angle: f64) -> TriMesh {
    let length = cyl_r * 2.0 + 4.0;
    let mut tun = make_cylinder(r_wire, length, 12);

    // Rotate so the cylinder lies along X
    rotate_y_mesh(&mut tun);

    // Centre it so it passes through the axis
    translate(&mut tun, -(length as f32 / 2.0), 0.0, z as f32);
    rotate_z(&mut tun, angle as f32);

    tun
}

/// Rotate all vertices 90 degrees about Y (Z -> X, X -> -Z).
fn rotate_y_mesh(mesh: &mut TriMesh) {
    for i in (0..mesh.positions.len()).step_by(3) {
        let x = mesh.positions[i];
        let z = mesh.positions[i + 2];
        mesh.positions[i] = z;
        mesh.positions[i + 2] = -x;
    }
}

/// Rotate all vertices about Y by an arbitrary angle.
fn rotate_y(mesh: &mut TriMesh, angle_rad: f32) {
    let (s, c) = (angle_rad.sin(), angle_rad.cos());
    for i in (0..mesh.positions.len()).step_by(3) {
        let x = mesh.positions[i];
        let z = mesh.positions[i + 2];
        mesh.positions[i] = x * c + z * s;
        mesh.positions[i + 2] = -x * s + z * c;
    }
}

// ─── STL export ────────────────────────────────────────────────────

/// Compute a face normal from three vertices (positions slice, indices).
fn face_normal(p: &[f32], i0: u32, i1: u32, i2: u32) -> [f32; 3] {
    let (i0, i1, i2) = (i0 as usize * 3, i1 as usize * 3, i2 as usize * 3);
    let ux = p[i1] - p[i0];
    let uy = p[i1 + 1] - p[i0 + 1];
    let uz = p[i1 + 2] - p[i0 + 2];
    let vx = p[i2] - p[i0];
    let vy = p[i2 + 1] - p[i0 + 1];
    let vz = p[i2 + 2] - p[i0 + 2];

    let nx = uy * vz - uz * vy;
    let ny = uz * vx - ux * vz;
    let nz = ux * vy - uy * vx;
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        [nx / len, ny / len, nz / len]
    }
}

/// Encode the mesh as a **binary STL** byte buffer.
pub fn to_binary_stl(mesh: &TriMesh) -> Vec<u8> {
    let tri_count = mesh.indices.len() / 3;
    // 80-byte header + 4-byte tri count + 50 bytes per triangle
    let size = 80 + 4 + tri_count * 50;
    let mut buf = Vec::with_capacity(size);

    // Header (80 bytes)
    let mut header = [0u8; 80];
    let label = b"Binary STL - W7HAK Coil Former";
    header[..label.len()].copy_from_slice(label);
    buf.extend_from_slice(&header);

    // Triangle count
    buf.extend_from_slice(&(tri_count as u32).to_le_bytes());

    for t in 0..tri_count {
        let i0 = mesh.indices[t * 3];
        let i1 = mesh.indices[t * 3 + 1];
        let i2 = mesh.indices[t * 3 + 2];

        let n = face_normal(&mesh.positions, i0, i1, i2);

        // normal
        buf.extend_from_slice(&n[0].to_le_bytes());
        buf.extend_from_slice(&n[1].to_le_bytes());
        buf.extend_from_slice(&n[2].to_le_bytes());

        // vertices
        for idx in [i0, i1, i2] {
            let base = idx as usize * 3;
            buf.extend_from_slice(&mesh.positions[base].to_le_bytes());
            buf.extend_from_slice(&mesh.positions[base + 1].to_le_bytes());
            buf.extend_from_slice(&mesh.positions[base + 2].to_le_bytes());
        }

        // attribute byte count
        buf.extend_from_slice(&0u16.to_le_bytes());
    }

    buf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_params_derived() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);

        assert!((d.rib_diam - 22.3).abs() < 1e-6);
        assert!((d.cylinder_diam - 19.1).abs() < 1e-6);
        assert!(d.calc_turns > 10.0);
        assert!(d.total_height > 0.0);
    }

    #[test]
    fn test_cylinder_mesh() {
        let cyl = make_cylinder(5.0, 10.0, 16);
        assert!(cyl.vertex_count() > 0);
        assert!(!cyl.indices.is_empty());
    }

    #[test]
    fn test_build_produces_mesh() {
        let p = CoilParams::default();
        let (mesh, _) = build_coil_former(&p);
        assert!(mesh.vertex_count() > 100);
        assert!(mesh.indices.len() > 100);
    }

    #[test]
    fn test_stl_export_size() {
        let p = CoilParams::default();
        let (mesh, _) = build_coil_former(&p);
        let stl = to_binary_stl(&mesh);
        let tri_count = mesh.indices.len() / 3;
        assert_eq!(stl.len(), 80 + 4 + tri_count * 50);
    }
}
