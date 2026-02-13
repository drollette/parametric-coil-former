//! Core coil-former geometry engine.
//!
//! Translates the parametric OpenSCAD / CadQuery model into pure Rust.
//! All units are millimetres.
//!
//! The body is built as a single integrated mesh: a hollow cylinder whose
//! outer radius varies to form the friction-rib profile and the helical
//! V-groove.  This avoids the CSG boolean operations the OpenSCAD version
//! relies on and produces a clean manifold mesh for both the 3-D preview
//! and the downloadable STL.

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
    #[allow(dead_code)]
    pub r_wire: f64,
    pub calc_turns: f64,
    pub winding_height: f64,
    pub total_height: f64,
    pub start_z: f64,
    pub end_z: f64,
    pub chamfer_h: f64,
    pub v_depth: f64,
    pub pitch: f64,
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
            calc_turns,
            winding_height,
            total_height,
            start_z,
            end_z,
            chamfer_h,
            v_depth,
            pitch: p.pitch,
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
    #[allow(dead_code)]
    pub fn append(&mut self, other: &TriMesh) {
        let base = (self.positions.len() / 3) as u32;
        self.positions.extend_from_slice(&other.positions);
        self.indices.extend(other.indices.iter().map(|i| i + base));
    }

    /// Number of vertices.
    #[cfg(test)]
    pub fn vertex_count(&self) -> usize {
        self.positions.len() / 3
    }
}

// ─── Cylinder generation (kept for tests / utility) ────────────────

/// Closed cylinder along the Z axis centred at origin.
#[allow(dead_code)]
pub fn make_cylinder(radius: f64, height: f64, segments: u32) -> TriMesh {
    let mut mesh = TriMesh::new();
    let segs = segments as usize;

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

        // side quad – outward-facing normals
        mesh.indices.extend_from_slice(&[b0, b1, t0]);
        mesh.indices.extend_from_slice(&[b1, t1, t0]);

        // bottom cap
        mesh.indices.extend_from_slice(&[bot_center, b0, b1]);

        // top cap
        mesh.indices.extend_from_slice(&[top_center, t1, t0]);
    }

    mesh
}

// ─── Transform helpers ─────────────────────────────────────────────

#[allow(dead_code)]
pub fn translate(mesh: &mut TriMesh, dx: f32, dy: f32, dz: f32) {
    for i in (0..mesh.positions.len()).step_by(3) {
        mesh.positions[i] += dx;
        mesh.positions[i + 1] += dy;
        mesh.positions[i + 2] += dz;
    }
}

#[allow(dead_code)]
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

/// Build the complete coil former mesh from parameters.
///
/// Produces a single integrated hollow cylinder whose outer surface
/// incorporates the friction-rib profile and helical V-groove.
pub fn build_coil_former(params: &CoilParams) -> (TriMesh, CoilDerived) {
    let d = CoilDerived::from_params(params);
    let mesh = build_body(&d);
    (mesh, d)
}

// ─── Integrated body mesh ──────────────────────────────────────────

/// Compute the set of z-levels used to tessellate the body.
fn compute_z_levels(d: &CoilDerived) -> Vec<f64> {
    // Base resolution: ~4 levels per mm
    let n_base = ((d.total_height * 4.0).ceil() as usize).max(200);

    let mut levels: Vec<f64> = (0..=n_base)
        .map(|i| d.total_height * (i as f64) / (n_base as f64))
        .collect();

    // Key transition points for crisp rib edges
    let bot_start = 2.0;
    let bot_chamfer_end = bot_start + d.chamfer_h;
    let bot_flat_end = bot_chamfer_end + 2.0;
    let bot_end = bot_flat_end + d.chamfer_h;

    let top_start = d.total_height - 8.0;
    let top_chamfer_end = top_start + d.chamfer_h;
    let top_flat_end = top_chamfer_end + 2.0;
    let top_end = top_flat_end + d.chamfer_h;

    for z in [
        bot_start,
        bot_chamfer_end,
        bot_flat_end,
        bot_end,
        top_start,
        top_chamfer_end,
        top_flat_end,
        top_end,
        d.start_z,
        d.end_z,
    ] {
        if z >= 0.0 && z <= d.total_height {
            levels.push(z);
        }
    }

    levels.sort_by(|a, b| a.partial_cmp(b).unwrap());
    levels.dedup_by(|a, b| (*a - *b).abs() < 0.01);
    levels
}

/// Outer radius at height `z` (rib profile, no groove).
fn rib_radius_at(d: &CoilDerived, z: f64) -> f64 {
    let cyl_r = d.cylinder_r;
    let rib_r = d.rib_diam / 2.0;

    // Bottom rib
    let bs = 2.0;
    let bc = bs + d.chamfer_h;
    let bf = bc + 2.0;
    let be = bf + d.chamfer_h;

    if z >= bs && z <= be {
        return if z <= bc {
            lerp(cyl_r, rib_r, (z - bs) / d.chamfer_h)
        } else if z <= bf {
            rib_r
        } else {
            lerp(rib_r, cyl_r, (z - bf) / d.chamfer_h)
        };
    }

    // Top rib
    let ts = d.total_height - 8.0;
    let tc = ts + d.chamfer_h;
    let tf = tc + 2.0;
    let te = tf + d.chamfer_h;

    if z >= ts && z <= te {
        return if z <= tc {
            lerp(cyl_r, rib_r, (z - ts) / d.chamfer_h)
        } else if z <= tf {
            rib_r
        } else {
            lerp(rib_r, cyl_r, (z - tf) / d.chamfer_h)
        };
    }

    cyl_r
}

/// V-groove depth at a surface point `(z, theta)`.
///
/// Uses the perpendicular distance from the point to the nearest helix
/// pass to determine whether the point lies inside the groove.
fn groove_depth_at(d: &CoilDerived, z: f64, theta: f64) -> f64 {
    if z < d.start_z || z > d.end_z {
        return 0.0;
    }

    // On the unwrapped cylinder surface the helix lines are straight and
    // spaced one pitch apart.  At angle θ the nearest groove centre is at
    // z = start_z + (−θ/(2π) + n)·pitch for some integer n.
    // Rearranging: (z − start_z + θ·pitch/(2π)) / pitch = n.
    let q = (z - d.start_z) / d.pitch + theta / (2.0 * PI);
    let z_mod = q.rem_euclid(1.0); // fractional part in [0, 1)
    let dz_frac = if z_mod <= 0.5 { z_mod } else { 1.0 - z_mod };
    let dz = dz_frac * d.pitch; // z-distance to nearest groove centre

    // Correct for helix angle to get true perpendicular distance
    let circ = 2.0 * PI * d.cylinder_r;
    let cos_alpha = circ / (circ * circ + d.pitch * d.pitch).sqrt();
    let d_perp = dz * cos_alpha;

    // 90° V-groove profile: linear ramp, depth = v_depth at centre
    if d_perp < d.v_depth {
        d.v_depth - d_perp
    } else {
        0.0
    }
}

/// Combined outer radius at `(z, theta)` including rib profile and groove.
fn radius_at(d: &CoilDerived, z: f64, theta: f64) -> f64 {
    let base_r = rib_radius_at(d, z);
    let groove = groove_depth_at(d, z, theta);
    (base_r - groove).max(d.center_bore_r + 0.1)
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t.clamp(0.0, 1.0) * (b - a)
}

/// Build the complete body as a single hollow mesh.
///
/// Outer surface: cylinder with rib profile and V-groove displacement.
/// Inner surface: constant centre-bore radius.
/// End caps: annular rings connecting outer ↔ inner.
fn build_body(d: &CoilDerived) -> TriMesh {
    let n_around = SEGMENTS as usize;
    let z_levels = compute_z_levels(d);
    let n_z = z_levels.len();
    let inner_r = d.center_bore_r;

    let mut mesh = TriMesh::new();

    // ── Outer surface vertices: [0 .. n_z * n_around) ──
    for &z in &z_levels {
        for i in 0..n_around {
            let theta = 2.0 * PI * (i as f64) / (n_around as f64);
            let r = radius_at(d, z, theta);
            mesh.positions.push((r * theta.cos()) as f32);
            mesh.positions.push((r * theta.sin()) as f32);
            mesh.positions.push(z as f32);
        }
    }

    // ── Inner surface vertices: [n_z*n_around .. 2*n_z*n_around) ──
    let inner_offset = n_z * n_around;
    for &z in &z_levels {
        for i in 0..n_around {
            let theta = 2.0 * PI * (i as f64) / (n_around as f64);
            mesh.positions.push((inner_r * theta.cos()) as f32);
            mesh.positions.push((inner_r * theta.sin()) as f32);
            mesh.positions.push(z as f32);
        }
    }

    // ── Outer side triangles (outward-facing normals) ──
    for j in 0..(n_z - 1) {
        for i in 0..n_around {
            let ni = (i + 1) % n_around;
            let v00 = (j * n_around + i) as u32;
            let v01 = (j * n_around + ni) as u32;
            let v10 = ((j + 1) * n_around + i) as u32;
            let v11 = ((j + 1) * n_around + ni) as u32;
            mesh.indices.extend_from_slice(&[v00, v01, v10]);
            mesh.indices.extend_from_slice(&[v01, v11, v10]);
        }
    }

    // ── Inner side triangles (inward-facing normals) ──
    for j in 0..(n_z - 1) {
        for i in 0..n_around {
            let ni = (i + 1) % n_around;
            let v00 = (inner_offset + j * n_around + i) as u32;
            let v01 = (inner_offset + j * n_around + ni) as u32;
            let v10 = (inner_offset + (j + 1) * n_around + i) as u32;
            let v11 = (inner_offset + (j + 1) * n_around + ni) as u32;
            mesh.indices.extend_from_slice(&[v00, v10, v01]);
            mesh.indices.extend_from_slice(&[v01, v10, v11]);
        }
    }

    // ── Bottom annular cap (normal pointing down, −z) ──
    for i in 0..n_around {
        let ni = (i + 1) % n_around;
        let o0 = i as u32;
        let o1 = ni as u32;
        let i0 = (inner_offset + i) as u32;
        let i1 = (inner_offset + ni) as u32;
        mesh.indices.extend_from_slice(&[o0, i0, o1]);
        mesh.indices.extend_from_slice(&[o1, i0, i1]);
    }

    // ── Top annular cap (normal pointing up, +z) ──
    let outer_top = (n_z - 1) * n_around;
    let inner_top = inner_offset + (n_z - 1) * n_around;
    for i in 0..n_around {
        let ni = (i + 1) % n_around;
        let o0 = (outer_top + i) as u32;
        let o1 = (outer_top + ni) as u32;
        let i0 = (inner_top + i) as u32;
        let i1 = (inner_top + ni) as u32;
        mesh.indices.extend_from_slice(&[o0, o1, i0]);
        mesh.indices.extend_from_slice(&[o1, i1, i0]);
    }

    mesh
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

    #[test]
    fn test_groove_depth_zero_outside_winding() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        assert_eq!(groove_depth_at(&d, 0.0, 0.0), 0.0);
        assert_eq!(groove_depth_at(&d, d.total_height, 0.0), 0.0);
    }

    #[test]
    fn test_groove_depth_max_on_helix() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        // On the helix centreline at z = start_z, theta = 0
        let depth = groove_depth_at(&d, d.start_z, 0.0);
        assert!((depth - d.v_depth).abs() < 0.01);
    }

    #[test]
    fn test_rib_radius_at_cylinder() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        // Outside rib zones → cylinder_r
        assert!((rib_radius_at(&d, 0.0) - d.cylinder_r).abs() < 1e-6);
        assert!((rib_radius_at(&d, d.total_height) - d.cylinder_r).abs() < 1e-6);
    }

    #[test]
    fn test_rib_radius_at_peak() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        let rib_r = d.rib_diam / 2.0;
        // Middle of bottom rib flat band
        let z_mid = 2.0 + d.chamfer_h + 1.0;
        assert!((rib_radius_at(&d, z_mid) - rib_r).abs() < 1e-6);
    }
}
