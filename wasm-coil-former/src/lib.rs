//! WASM Coil Former – browser-facing entry point.
//!
//! Compiled to WebAssembly via `wasm-pack build --target web`.
//! Exposes functions for the JS frontend to:
//!   1. Generate a 3D preview mesh (positions + indices)
//!   2. Export a binary STL file

mod geometry;

use geometry::{build_coil_former, to_binary_stl, CoilDerived, CoilParams};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

/// Version string baked in at compile time.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// JSON-serialisable result returned to JavaScript.
#[derive(Serialize, Deserialize)]
pub struct CoilResult {
    /// Flat f32 vertex positions [x,y,z, x,y,z, …]
    pub positions: Vec<f32>,
    /// Triangle indices
    pub indices: Vec<u32>,
    /// Computed turns
    pub turns: f64,
    /// Winding height in mm
    pub coil_height: f64,
    /// Total former height in mm
    pub total_height: f64,
    /// Cylinder outer diameter
    pub cylinder_diam: f64,
}

/// Generate the coil former mesh and return the result as a JSON string.
///
/// Called from JavaScript with the five user-facing parameters.
#[wasm_bindgen]
pub fn generate_coil(
    wire_len: f64,
    wire_diam: f64,
    pvc_inner_diam: f64,
    pitch: f64,
    rib_clearance: f64,
) -> String {
    let params = CoilParams {
        wire_len,
        wire_diam,
        pvc_inner_diam,
        pitch,
        rib_clearance,
    };

    let (mesh, derived) = build_coil_former(&params);

    let result = CoilResult {
        positions: mesh.positions,
        indices: mesh.indices,
        turns: derived.calc_turns,
        coil_height: derived.winding_height,
        total_height: derived.total_height,
        cylinder_diam: derived.cylinder_diam,
    };

    serde_json::to_string(&result).unwrap_or_else(|e| format!("{{\"error\":\"{}\"}}", e))
}

/// Generate a binary STL file and return it as a `Uint8Array` for download.
#[wasm_bindgen]
pub fn export_stl(
    wire_len: f64,
    wire_diam: f64,
    pvc_inner_diam: f64,
    pitch: f64,
    rib_clearance: f64,
) -> Vec<u8> {
    let params = CoilParams {
        wire_len,
        wire_diam,
        pvc_inner_diam,
        pitch,
        rib_clearance,
    };

    let (mesh, _) = build_coil_former(&params);
    to_binary_stl(&mesh)
}

/// Return only the computed dimensions as JSON (no mesh data).
/// Useful for a lightweight info panel update without re-rendering.
#[wasm_bindgen]
pub fn compute_info(
    wire_len: f64,
    wire_diam: f64,
    pvc_inner_diam: f64,
    pitch: f64,
    rib_clearance: f64,
) -> String {
    let params = CoilParams {
        wire_len,
        wire_diam,
        pvc_inner_diam,
        pitch,
        rib_clearance,
    };
    let d = CoilDerived::from_params(&params);

    serde_json::json!({
        "turns": d.calc_turns,
        "coil_height": d.winding_height,
        "total_height": d.total_height,
        "cylinder_diam": d.cylinder_diam,
        "rib_diam": d.rib_diam,
        "center_bore_diam": d.center_bore_r * 2.0,
        "v_depth": d.v_depth,
    })
    .to_string()
}
