# Spectral BRDF Simulation and Color Reconstruction 

This project computes the reflected color of a material when illuminated by a user-defined light source using spectral albedo data provided as CSV files.

## Overview

The model applies:

* CIE 1931 spectral integration
* Simplified Phong BRDF
* Delaunay interpolation in xy chromaticity space
* XYZ → xyY → sRGB conversions

Outputs include both numeric XYZ values and a visual color patch.
Although the context is optical, the core of the project is highly quantitative: it involves numerical integration, coordinate transforms, interpolation on irregular meshes, and geometric modeling.

## Usage

1. **Select a CSV albedo file** from the provided folder, e.g.:

```
led_4000K_white.csv
```

2. **Call the main function**:

```python
couleur_reflexion(
    albedo_csv,           # path to CSV file
    couleur_source_XYZ,   # XYZ vector of light source
    angle_source,         # angle of source in degrees
    angle_capteur,        # angle of captor in degrees
    n,                    # Phong specular exponent
    ks                    # specular coefficient
)
```

### Example:

```python
couleur_reflexion(
    "led_4000K_white.csv",
    [1.05, 1.00, 1.20],  # light source XYZ inside mesh
    0,                    # source angle
    10,                   # captor angle
    20,                   # Phong exponent
    0.05                  # specular coefficient
)
```

The function will:

1. Convert albedo spectrum to XYZ
2. Apply BRDF correction
3. Interpolate chromaticity in Delaunay mesh
4. Convert to XYZ and sRGB
5. Display the reflected color patch

## Scientific & Quantitative Features

* Spectral integration using CIE 1931 2° Standard Observer
* Phong BRDF modeling with user-defined parameters
* Delaunay triangulation and barycentric interpolation in xy space
* Full XYZ ↔ xyY ↔ sRGB pipeline with gamma correction
* Reproducible, vectorized computation suitable for quantitative modeling

## Output

* Displayed color patch representing reflection
* XYZ numeric values shown in plot title

Simply pick a CSV file and run the function with the parameters you want.
