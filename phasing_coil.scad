// UNIVERSAL Phasing Coil Former - FULL MASTER SCRIPT
// Inputs: Wire Length, Wire Diam, PVC Inner Diam, Pitch

// ==========================================
// --- MAIN USER INPUTS (CONTROL PANEL) ---
// ==========================================                                           --- ORIGINAL VALUES
wire_len = 668;           // 1. Target wire length (Electrical phase length)            --- 668
wire_diam = 3.2;          // 2. Wire diameter (e.g., 3.2 for RG58 core, 1.6 for bare)   --- 3.2
pvc_inner_diam = 23.3;    // 3. Exact inside diameter of your PVC pipe                  --- 23.3            
pitch = 8.9;              // 4. Vertical distance between each wrap                     --- 8.9

// --- SECONDARY SETTINGS ---
rib_clearance = 8;        // Solid plastic space at top and bottom for the ribs
$fn = 60;                 // 3D rendering smoothness

// --- PERFORMANCE TUNING ---
groove_steps_per_turn = 24;  // Reduce for faster render (24-50), increase for smoother groove
cylinder_faces = 72;         // Faces for main cylinders (72-100)

// ==========================================
// --- AUTOMATIC CALCULATIONS (Do not touch) ---
// ==========================================

// Auto-size the ribs and main body to ensure wire sits perfectly flush
rib_diam = pvc_inner_diam; 
cylinder_diam = rib_diam - wire_diam;     

// Auto-size the center hollow shaft to fit the wire exactly
center_bore_diam = wire_diam + 0.2;       

// Helical Math
circumference = 3.14159 * cylinder_diam;
calc_turns = sqrt(pow(wire_len, 2) / (pow(circumference, 2) + pow(pitch, 2)));
wh = calc_turns * pitch;
total_height = wh + (rib_clearance * 2);

start_z = rib_clearance;
end_z = start_z + wh;
chamfer_h = (rib_diam - cylinder_diam) / 2; 

// The exact depth a 90-deg V must be to cradle 50% of the wire
v_tip_depth = (wire_diam / 2) * sqrt(2); 

// Print out the vital stats to the console
echo("==============================");
echo(str("Wire Length: ", wire_len, " mm"));
echo(str("Wire Diam: ", wire_diam, " mm"));
echo(str("PVC Inner Diam: ", pvc_inner_diam, " mm"));
echo(str("Pitch: ", pitch, " mm"));
echo("---");
echo(str("Calculated Turns: ", calc_turns)); 
echo(str("TOTAL FORMER HEIGHT: ", total_height, " mm"));
echo("==============================");

difference() {
    union() {
        // Main Body
        cylinder(h = total_height, d = cylinder_diam, $fn=cylinder_faces);
        
        // --- Bottom Rib ---
        translate([0, 0, 2]) {
            cylinder(h = chamfer_h, d1 = cylinder_diam, d2 = rib_diam, $fn=cylinder_faces);
            translate([0, 0, chamfer_h]) cylinder(h = 2, d = rib_diam, $fn=cylinder_faces);
            translate([0, 0, chamfer_h + 2]) cylinder(h = chamfer_h, d1 = rib_diam, d2 = cylinder_diam, $fn=cylinder_faces);
        }

        // --- Top Rib (Auto-positions based on calculated height) ---
        translate([0, 0, total_height - 8]) {
            cylinder(h = chamfer_h, d1 = cylinder_diam, d2 = rib_diam, $fn=cylinder_faces);
            translate([0, 0, chamfer_h]) cylinder(h = 2, d = rib_diam, $fn=cylinder_faces);
            translate([0, 0, chamfer_h + 2]) cylinder(h = chamfer_h, d1 = rib_diam, d2 = cylinder_diam, $fn=cylinder_faces);
        }
    }

    // --- 1. Center Void ---
    translate([0, 0, -1])
        cylinder(h = total_height + 2, d = center_bore_diam, $fn=cylinder_faces);

    // --- 2. Flawless 3D V-Groove ---
    translate([0, 0, start_z])
        true_3d_v_groove(h=wh, turns=calc_turns, d=cylinder_diam, depth=v_tip_depth);

    // --- 3. Swept Entry/Exit Tunnels ---
    hull() {
        translate([0, 0, start_z - 4]) sphere(d=wire_diam);
        translate([cylinder_diam/2, 0, start_z]) sphere(d=wire_diam);
        rotate([0, 0, -10]) 
            translate([cylinder_diam/2, 0, start_z + wh*(10/(360*calc_turns))]) 
            sphere(d=wire_diam);
    }

    hull() {
        translate([0, 0, end_z + 4]) sphere(d=wire_diam);
        rotate([0, 0, -360 * calc_turns]) 
            translate([cylinder_diam/2, 0, end_z]) sphere(d=wire_diam);
        rotate([0, 0, -(360 * calc_turns - 10)]) 
            translate([cylinder_diam/2, 0, end_z - wh*(10/(360*calc_turns))]) 
            sphere(d=wire_diam);
    }
}

// --- True 3D Sweep Engine ---
module true_3d_v_groove(h, turns, d, depth) {
    total_steps = round(turns * groove_steps_per_turn);
    dz = h / total_steps;
    da = -360 * turns / total_steps;

    for (i = [0 : total_steps - 1]) {
        hull() {
            cutter_pos(i * dz, i * da, d, depth);
            cutter_pos((i+1) * dz, (i+1) * da, d, depth);
        }
    }
}

module cutter_pos(z, a, d, depth) {
    translate([0, 0, z])
    rotate([0, 0, a])
    translate([d/2 - depth, 0, 0]) 
    rotate([0, 90, 0])
    cylinder(r1=0, r2=depth + 1, h=depth + 1, $fn=16);
}