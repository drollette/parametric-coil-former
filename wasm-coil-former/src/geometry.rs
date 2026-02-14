//! Core coil-former geometry engine.
//!
//! Translates the parametric OpenSCAD / CadQuery model into pure Rust.
//! All units are millimetres.
//!
//! The body is built as a single integrated mesh: a hollow cylinder whose
//! outer radius varies to form the friction-rib profile, helical V-groove,
//! and wire entry/exit holes with chamfered edges.

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
    /// Centre bore diameter in mm.  Pass 0 (or negative) for automatic
    /// sizing (wire_diam + 0.2mm clearance).
    pub center_bore_diam: f64,
}

impl Default for CoilParams {
    fn default() -> Self {
        Self {
            wire_len: 668.0,
            wire_diam: 3.2,
            pvc_inner_diam: 22.3,
            pitch: 8.9,
            rib_clearance: 8.0,
            center_bore_diam: 0.0, // auto
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
    pub calc_turns: f64,
    pub winding_height: f64,
    pub total_height: f64,
    pub start_z: f64,
    pub end_z: f64,
    pub chamfer_h: f64,
    pub v_depth: f64,
    pub pitch: f64,
    /// Entry-hole angle (always 0).
    pub entry_angle: f64,
    /// Exit-hole angle, normalised to [0, 2π).
    #[allow(dead_code)]
    pub exit_angle: f64,
    /// Channel extension distance from groove to bore (wall thickness + margin).
    pub channel_extension: f64,
}

impl CoilDerived {
    pub fn from_params(p: &CoilParams) -> Self {
        let rib_diam = p.pvc_inner_diam;
        let cylinder_diam = rib_diam - p.wire_diam;
        let cylinder_r = cylinder_diam / 2.0;
        let center_bore_diam = if p.center_bore_diam > 0.0 {
            p.center_bore_diam
        } else {
            p.wire_diam + 0.2
        };
        let center_bore_r = center_bore_diam / 2.0;
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

        let entry_angle = 0.0;
        let exit_angle = ((-2.0 * PI * calc_turns) % (2.0 * PI) + 2.0 * PI) % (2.0 * PI);

        // Channel extension: axial distance the wire slot extends beyond the
        // winding zone.  Half the bore diameter gives a proportional opening
        // that won't cut through ribs or end caps.
        let channel_extension = center_bore_r;

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
            entry_angle,
            exit_angle,
            channel_extension,
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

        mesh.indices.extend_from_slice(&[b0, b1, t0]);
        mesh.indices.extend_from_slice(&[b1, t1, t0]);
        mesh.indices.extend_from_slice(&[bot_center, b0, b1]);
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

const SEGMENTS: u32 = 128;

/// Build the complete coil former mesh from parameters.
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

    // Channel boundary z-levels for clean tunnel edges
    let z_fillet = d.r_wire;
    let entry_z_bore = d.start_z - d.channel_extension;
    let exit_z_bore = d.end_z + d.channel_extension;

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
        // Entry channel boundaries
        entry_z_bore - z_fillet,
        entry_z_bore,
        d.start_z + z_fillet,
        // Exit channel boundaries
        d.end_z - z_fillet,
        exit_z_bore,
        exit_z_bore + z_fillet,
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

// ─── Continuous wire path geometry ────────────────────────────────

/// A point along the continuous wire path with position and orientation.
#[derive(Debug, Clone, Copy)]
struct PathPoint {
    /// Position in 3D space (x, y, z)
    position: [f64; 3],
    /// Unit tangent vector (direction of wire travel)
    #[allow(dead_code)]
    tangent: [f64; 3],
    /// Distance parameter along path (0 = entry, total_path_length = exit)
    #[allow(dead_code)]
    distance: f64,
}

/// Cubic Bézier curve for smooth wire path transitions.
/// Creates a 90° elbow from vertical bore to horizontal groove surface.
fn cubic_bezier(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], t: f64) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let s = 1.0 - t;
    let s2 = s * s;
    let s3 = s2 * s;
    let t2 = t * t;
    let t3 = t2 * t;

    [
        s3 * p0[0] + 3.0 * s2 * t * p1[0] + 3.0 * s * t2 * p2[0] + t3 * p3[0],
        s3 * p0[1] + 3.0 * s2 * t * p1[1] + 3.0 * s * t2 * p2[1] + t3 * p3[1],
        s3 * p0[2] + 3.0 * s2 * t * p1[2] + 3.0 * s * t2 * p2[2] + t3 * p3[2],
    ]
}

/// Tangent (derivative) of cubic Bézier curve.
fn cubic_bezier_tangent(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    t: f64,
) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let s = 1.0 - t;

    let dx = 3.0 * s * s * (p1[0] - p0[0])
        + 6.0 * s * t * (p2[0] - p1[0])
        + 3.0 * t * t * (p3[0] - p2[0]);
    let dy = 3.0 * s * s * (p1[1] - p0[1])
        + 6.0 * s * t * (p2[1] - p1[1])
        + 3.0 * t * t * (p3[1] - p2[1]);
    let dz = 3.0 * s * s * (p1[2] - p0[2])
        + 6.0 * s * t * (p2[2] - p1[2])
        + 3.0 * t * t * (p3[2] - p2[2]);

    normalize([dx, dy, dz])
}

/// Normalize a 3D vector to unit length.
fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Dot product of two 3D vectors.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the helix tangent vector at a given turn fraction `t` ∈ [0, 1].
///
/// The helix is parameterised as:
///   θ(t) = entry_angle − 2π·turns·t
///   pos   = (R cos θ, R sin θ, start_z + t·winding_height)
///
/// Returns the *unnormalised* derivative so the caller can normalise or
/// scale it as needed for Bézier handle length.
fn helix_tangent(d: &CoilDerived, t: f64) -> [f64; 3] {
    let theta = d.entry_angle - 2.0 * PI * d.calc_turns * t;
    let omega = -2.0 * PI * d.calc_turns; // dθ/dt
    [
        -d.cylinder_r * theta.sin() * omega,
        d.cylinder_r * theta.cos() * omega,
        d.winding_height,
    ]
}

/// Helix position at turn fraction `t` ∈ [0, 1].
fn helix_pos(d: &CoilDerived, t: f64) -> [f64; 3] {
    let theta = d.entry_angle - 2.0 * PI * d.calc_turns * t;
    [
        d.cylinder_r * theta.cos(),
        d.cylinder_r * theta.sin(),
        d.start_z + t * d.winding_height,
    ]
}

/// Scale a vector by a scalar.
fn vec_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Add two vectors.
fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two vectors (a - b).
fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Length of a vector.
fn vec_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Generate the complete continuous wire path from entry to exit.
///
/// The path consists of 5 segments:
/// 1. Vertical bore (bottom → transition start)
/// 2. Bézier elbow  (bore → helix start, tangent-matched to helix)
/// 3. Helical winding
/// 4. Bézier elbow  (helix end → bore, tangent-matched to helix)
/// 5. Vertical bore  (transition end → top)
///
/// Key fixes vs. the previous version:
/// - Bézier endpoints sit 0.5 mm *outside* the former surface so the
///   subtraction cutter overshoots and leaves no ghost faces.
/// - Control points are derived from the actual helix tangent at t=0 / t=1
///   so the groove and hole meet seamlessly (no angle mismatch).
fn generate_wire_path(d: &CoilDerived) -> Vec<PathPoint> {
    let mut path = Vec::new();

    // ── Overcut margin: extend cutter 0.5 mm past the outer surface ──
    let overcut = 0.5_f64;

    // Helix tangent at start (t = 0) and end (t = 1)
    let helix_tan_start = normalize(helix_tangent(d, 0.0)); // direction wire exits bore
    let helix_tan_end = normalize(helix_tangent(d, 1.0)); // direction wire enters bore

    // Helix positions at start / end, pushed 0.5 mm outward for overcut
    let helix_start = helix_pos(d, 0.0);
    let helix_end = helix_pos(d, 1.0);

    // Radial outward direction at helix start / end (unit)
    let radial_start = normalize([helix_start[0], helix_start[1], 0.0]);
    let radial_end = normalize([helix_end[0], helix_end[1], 0.0]);

    // Overcut surface points (slightly outside the former)
    let surface_start = vec_add(helix_start, vec_scale(radial_start, overcut));
    let surface_end = vec_add(helix_end, vec_scale(radial_end, overcut));

    // Handle length for Bézier curves – controls how "wide" the elbow is.
    // Use the wall thickness (cylinder_r − bore_r) so the arc clears
    // the wall comfortably.
    let wall = d.cylinder_r - d.center_bore_r;
    let handle_len = wall * 0.6;

    // ──── SEGMENT 1: Vertical bore (entry) ────
    // Start below the winding zone so the slot exits the bottom face.
    let bore_entry_z = d.start_z - d.channel_extension;
    // The bore-side end of the transition (where curve begins)
    let bore_entry_top_z = helix_start[2];
    let bore_entry_top = [0.0, 0.0, bore_entry_top_z];

    let seg1_len = (bore_entry_top_z - bore_entry_z).abs();
    let n_bore_entry = ((seg1_len / 0.3).ceil() as usize).max(2);
    for i in 0..=n_bore_entry {
        let t = i as f64 / n_bore_entry as f64;
        let z = bore_entry_z + t * (bore_entry_top_z - bore_entry_z);
        path.push(PathPoint {
            position: [0.0, 0.0, z],
            tangent: [0.0, 0.0, 1.0],
            distance: t * seg1_len,
        });
    }
    let mut dist_so_far = seg1_len;

    // ──── SEGMENT 2: Bézier elbow (bore → helix start) ────
    // P0 = bore top,  tangent = +Z
    // P3 = overcut helix start,  tangent = helix_tan_start
    //
    // P1 = P0 + handle_len * (0,0,1)        — continues upward from bore
    // P2 = P3 − handle_len * helix_tan_start — approaches helix smoothly
    let p0 = bore_entry_top;
    let p3 = surface_start;
    let p1 = vec_add(p0, [0.0, 0.0, handle_len]);
    let p2 = vec_sub(p3, vec_scale(helix_tan_start, handle_len));

    let elbow_len = vec_len(vec_sub(p3, p0)) * 1.2; // rough arc length
    let n_elbow = ((elbow_len / 0.3).ceil() as usize).max(8);
    for i in 1..=n_elbow {
        let t = i as f64 / n_elbow as f64;
        path.push(PathPoint {
            position: cubic_bezier(p0, p1, p2, p3, t),
            tangent: cubic_bezier_tangent(p0, p1, p2, p3, t),
            distance: dist_so_far + t * elbow_len,
        });
    }
    dist_so_far += elbow_len;

    // ──── SEGMENT 3: Helical winding ────
    // Sample densely (every ~0.3 mm of arc) for a smooth groove.
    let circ = 2.0 * PI * d.cylinder_r;
    let helix_arc = (circ * circ + d.pitch * d.pitch).sqrt() * d.calc_turns;
    let n_helix = ((helix_arc / 0.3).ceil() as usize).max(64);
    for i in 1..=n_helix {
        let t = i as f64 / n_helix as f64;
        let pos = helix_pos(d, t);
        let tan = normalize(helix_tangent(d, t));
        // Push the path point outward by the overcut so the cutter
        // extends through the outer wall.
        let radial = normalize([pos[0], pos[1], 0.0]);
        path.push(PathPoint {
            position: vec_add(pos, vec_scale(radial, overcut)),
            tangent: tan,
            distance: dist_so_far + t * helix_arc,
        });
    }
    dist_so_far += helix_arc;

    // ──── SEGMENT 4: Bézier elbow (helix end → bore) ────
    // P0 = overcut helix end,  tangent = helix_tan_end
    // P3 = bore point,         tangent = +Z (upward into bore)
    let bore_exit_bottom_z = helix_end[2];
    let bore_exit_bottom = [0.0, 0.0, bore_exit_bottom_z];

    let q0 = surface_end;
    let q3 = bore_exit_bottom;
    let q1 = vec_add(q0, vec_scale(helix_tan_end, handle_len));
    let q2 = vec_sub(q3, [0.0, 0.0, handle_len]); // approach bore from above

    let elbow_len2 = vec_len(vec_sub(q3, q0)) * 1.2;
    let n_elbow2 = ((elbow_len2 / 0.3).ceil() as usize).max(8);
    for i in 1..=n_elbow2 {
        let t = i as f64 / n_elbow2 as f64;
        path.push(PathPoint {
            position: cubic_bezier(q0, q1, q2, q3, t),
            tangent: cubic_bezier_tangent(q0, q1, q2, q3, t),
            distance: dist_so_far + t * elbow_len2,
        });
    }
    dist_so_far += elbow_len2;

    // ──── SEGMENT 5: Vertical bore (exit) ────
    let bore_exit_z = d.end_z + d.channel_extension;
    let seg5_len = (bore_exit_z - bore_exit_bottom_z).abs();
    let n_bore_exit = ((seg5_len / 0.3).ceil() as usize).max(2);
    for i in 1..=n_bore_exit {
        let t = i as f64 / n_bore_exit as f64;
        let z = bore_exit_bottom_z + t * (bore_exit_z - bore_exit_bottom_z);
        path.push(PathPoint {
            position: [0.0, 0.0, z],
            tangent: [0.0, 0.0, 1.0],
            distance: dist_so_far + t * seg5_len,
        });
    }

    path
}

/// Calculate signed distance from a point to the continuous wire path.
///
/// Returns negative if inside the wire volume, positive if outside.
/// Uses swept-sphere approach with slight overlap (1.05×) to ensure manifold geometry.
fn distance_to_wire_path(d: &CoilDerived, path: &[PathPoint], point: [f64; 3]) -> f64 {
    let wire_r = d.r_wire * 1.05; // 5% overlap margin for clean boolean

    let mut min_dist = f64::INFINITY;

    for window in path.windows(2) {
        let p0 = window[0].position;
        let p1 = window[1].position;

        // Vector from p0 to p1
        let edge = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let edge_len_sq = edge[0] * edge[0] + edge[1] * edge[1] + edge[2] * edge[2];

        if edge_len_sq < 1e-12 {
            continue;
        }

        // Vector from p0 to point
        let to_point = [point[0] - p0[0], point[1] - p0[1], point[2] - p0[2]];

        // Project point onto edge
        let t = (dot(to_point, edge) / edge_len_sq).clamp(0.0, 1.0);

        // Closest point on edge
        let closest = [
            p0[0] + t * edge[0],
            p0[1] + t * edge[1],
            p0[2] + t * edge[2],
        ];

        // Distance from point to closest point on edge
        let dx = point[0] - closest[0];
        let dy = point[1] - closest[1];
        let dz = point[2] - closest[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        min_dist = min_dist.min(dist);
    }

    min_dist - wire_r
}

// ─── Combined radius ───────────────────────────────────────────────

/// Combined outer radius at `(z, theta)` including rib profile and continuous
/// wire path subtraction (groove + transitions).
fn radius_at(d: &CoilDerived, z: f64, theta: f64, wire_path: &[PathPoint]) -> f64 {
    let base_r = rib_radius_at(d, z);

    // Convert cylindrical (r, theta, z) to cartesian for distance query
    let r_test = base_r;
    let x = r_test * theta.cos();
    let y = r_test * theta.sin();
    let point = [x, y, z];

    // Check if point is inside the continuous wire path
    let dist = distance_to_wire_path(d, wire_path, point);

    if dist < 0.0 {
        // Inside wire volume - carve down to bore or minimum safe distance
        // Binary search for the radius where we exit the wire volume
        let mut r_min = d.center_bore_r;
        let mut r_max = base_r;

        for _ in 0..10 {
            // 10 iterations gives ~0.1% precision
            let r_mid = (r_min + r_max) / 2.0;
            let x_mid = r_mid * theta.cos();
            let y_mid = r_mid * theta.sin();
            let dist_mid = distance_to_wire_path(d, wire_path, [x_mid, y_mid, z]);

            if dist_mid < 0.0 {
                r_max = r_mid; // Still inside, reduce radius
            } else {
                r_min = r_mid; // Outside, increase radius
            }
        }

        r_min.max(d.center_bore_r)
    } else {
        base_r.max(d.center_bore_r + 0.1)
    }
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t.clamp(0.0, 1.0) * (b - a)
}

/// Build the complete body as a single hollow mesh.
fn build_body(d: &CoilDerived) -> TriMesh {
    let n_around = SEGMENTS as usize;
    let z_levels = compute_z_levels(d);
    let n_z = z_levels.len();
    let inner_r = d.center_bore_r;

    // Generate continuous wire path once
    let wire_path = generate_wire_path(d);

    let mut mesh = TriMesh::new();

    // ── Outer surface vertices ──
    for &z in &z_levels {
        for i in 0..n_around {
            let theta = 2.0 * PI * (i as f64) / (n_around as f64);
            let r = radius_at(d, z, theta, &wire_path);
            mesh.positions.push((r * theta.cos()) as f32);
            mesh.positions.push((r * theta.sin()) as f32);
            mesh.positions.push(z as f32);
        }
    }

    // ── Inner surface vertices ──
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

    // ── Bottom annular cap (−z) ──
    for i in 0..n_around {
        let ni = (i + 1) % n_around;
        let o0 = i as u32;
        let o1 = ni as u32;
        let i0 = (inner_offset + i) as u32;
        let i1 = (inner_offset + ni) as u32;
        mesh.indices.extend_from_slice(&[o0, i0, o1]);
        mesh.indices.extend_from_slice(&[o1, i0, i1]);
    }

    // ── Top annular cap (+z) ──
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
    let size = 80 + 4 + tri_count * 50;
    let mut buf = Vec::with_capacity(size);

    let mut header = [0u8; 80];
    let label = b"Binary STL - W7HAK Coil Former";
    header[..label.len()].copy_from_slice(label);
    buf.extend_from_slice(&header);
    buf.extend_from_slice(&(tri_count as u32).to_le_bytes());

    for t in 0..tri_count {
        let i0 = mesh.indices[t * 3];
        let i1 = mesh.indices[t * 3 + 1];
        let i2 = mesh.indices[t * 3 + 2];

        let n = face_normal(&mesh.positions, i0, i1, i2);
        buf.extend_from_slice(&n[0].to_le_bytes());
        buf.extend_from_slice(&n[1].to_le_bytes());
        buf.extend_from_slice(&n[2].to_le_bytes());

        for idx in [i0, i1, i2] {
            let base = idx as usize * 3;
            buf.extend_from_slice(&mesh.positions[base].to_le_bytes());
            buf.extend_from_slice(&mesh.positions[base + 1].to_le_bytes());
            buf.extend_from_slice(&mesh.positions[base + 2].to_le_bytes());
        }
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
    fn test_custom_bore_diameter() {
        let p = CoilParams {
            center_bore_diam: 5.0,
            ..CoilParams::default()
        };
        let d = CoilDerived::from_params(&p);
        assert!((d.center_bore_r - 2.5).abs() < 1e-6);
    }

    #[test]
    fn test_auto_bore_diameter() {
        let p = CoilParams::default(); // center_bore_diam = 0 → auto
        let d = CoilDerived::from_params(&p);
        let expected = (3.2 + 0.2) / 2.0; // (wire_diam + 0.2) / 2
        assert!((d.center_bore_r - expected).abs() < 1e-6);
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
    fn test_wire_path_generation() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        let path = generate_wire_path(&d);
        // Path should have multiple points
        assert!(path.len() > 100, "Path has {} points", path.len());
        // First point should be near entry bore start
        let entry_bore_start_z = d.start_z - d.channel_extension;
        eprintln!(
            "First point z: {}, expected: {}",
            path[0].position[2], entry_bore_start_z
        );
        eprintln!(
            "Last point z: {}, expected: {}",
            path[path.len() - 1].position[2],
            d.end_z + d.channel_extension
        );
        assert!(
            (path[0].position[2] - entry_bore_start_z).abs() < 1.0,
            "First z={}, expected={}",
            path[0].position[2],
            entry_bore_start_z
        );
        // Last point should be near exit bore end
        let exit_bore_end_z = d.end_z + d.channel_extension;
        assert!(
            (path[path.len() - 1].position[2] - exit_bore_end_z).abs() < 1.0,
            "Last z={}, expected={}",
            path[path.len() - 1].position[2],
            exit_bore_end_z
        );
    }

    #[test]
    fn test_distance_to_wire_path() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        let path = generate_wire_path(&d);

        // Point at bore center should be inside wire path
        let bore_point = [0.0, 0.0, d.start_z];
        let dist = distance_to_wire_path(&d, &path, bore_point);
        assert!(dist < 0.0, "Bore center should be inside wire path");

        // Point far from path should be outside
        let far_point = [d.cylinder_r * 2.0, 0.0, d.total_height / 2.0];
        let dist_far = distance_to_wire_path(&d, &path, far_point);
        assert!(dist_far > 0.0, "Far point should be outside wire path");
    }

    #[test]
    fn test_rib_radius_at_cylinder() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        assert!((rib_radius_at(&d, 0.0) - d.cylinder_r).abs() < 1e-6);
        assert!((rib_radius_at(&d, d.total_height) - d.cylinder_r).abs() < 1e-6);
    }

    #[test]
    fn test_rib_radius_at_peak() {
        let p = CoilParams::default();
        let d = CoilDerived::from_params(&p);
        let rib_r = d.rib_diam / 2.0;
        let z_mid = 2.0 + d.chamfer_h + 1.0;
        assert!((rib_radius_at(&d, z_mid) - rib_r).abs() < 1e-6);
    }

    #[test]
    fn test_bezier_curve() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [0.0, 0.0, 1.0];
        let p2 = [1.0, 0.0, 1.0];
        let p3 = [1.0, 0.0, 0.0];

        // At t=0, should be at p0
        let start = cubic_bezier(p0, p1, p2, p3, 0.0);
        assert!((start[0] - p0[0]).abs() < 1e-6);
        assert!((start[1] - p0[1]).abs() < 1e-6);
        assert!((start[2] - p0[2]).abs() < 1e-6);

        // At t=1, should be at p3
        let end = cubic_bezier(p0, p1, p2, p3, 1.0);
        assert!((end[0] - p3[0]).abs() < 1e-6);
        assert!((end[1] - p3[1]).abs() < 1e-6);
        assert!((end[2] - p3[2]).abs() < 1e-6);
    }
}
