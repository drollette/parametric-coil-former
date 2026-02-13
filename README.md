# Phasing Coil Former Generator

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/drollette/parametric-coil-former/master)

An OpenSCAD parametric generator for 3D-printable phasing coil formers, designed for HAM radio antenna phasing lines.

## What It Does

This script generates a cylindrical coil former with a precision V-groove helix that guides wire placement. The former is sized to friction-fit inside standard PVC pipe, making it easy to create weatherproof phasing coils for antenna arrays.

## Parameters

Edit these values at the top of `phasing_coil.scad`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `wire_len` | Target wire length in mm (electrical phase length) | 668 |
| `wire_diam` | Wire diameter in mm (3.2 for RG58 core, 1.6 for bare wire) | 3.2 |
| `pvc_inner_diam` | Inside diameter of your PVC pipe in mm | 23.3 |
| `pitch` | Vertical spacing between wraps in mm | 8.9 |

## Features

- Automatic calculation of turns based on wire length and pitch
- V-groove sized to cradle 50% of wire diameter
- Swept entry/exit tunnels for clean wire routing
- Center bore for wire pass-through
- Top and bottom ribs for friction fit in PVC pipe

## Usage

1. Open `phasing_coil.scad` in OpenSCAD
2. Adjust the parameters for your wire and PVC pipe
3. Render (F6) and export as STL
4. Print with ~20% infill

## Try It Online

Click the Binder badge above or use the link below to launch an interactive Jupyter environment with CadQuery and jupyter-cadquery pre-installed:

https://mybinder.org/v2/gh/drollette/parametric-coil-former/master

Once the environment starts, open `viewer_check.ipynb` to verify that CadQuery rendering is working.

## License

MIT License - See [LICENSE](LICENSE) file

## Author

W7HAK
