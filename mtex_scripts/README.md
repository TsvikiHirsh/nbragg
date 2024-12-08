# nbragg Texture Analysis Scripts

## Overview
This folder contains MATLAB scripts for processing diffraction and pole figure data, specifically designed to work with the MTEX toolbox and prepare data for analysis with the nbragg CrossSection method.

## Prerequisites
- MATLAB
- MTEX Toolbox
- Crystal symmetry data for the material of interest

## Files Description

### Data Files
- `RPV_alpha.xpc`: Pole figure data for alpha phase (likely iron)
- `RPV_gamma.xpc`: Pole figure data for gamma phase (likely iron)

### MATLAB Scripts

#### Pole Figure Loading and Processing
- `loadPoleFigure_beartex.m`: Custom function to import pole figure data from BeaTex file format
  - Usage: `pf = loadPoleFigure_beartex(filename)`
  - Supports loading experimental pole figure data 

#### Orientation Component Extraction
- `alpha_component_extraction.m`: Extract texture components for alpha phase
  - Configures MTEX preferences
  - Loads alpha phase pole figure data
  - Prepares data for further analysis

- `gamma_component_extraction.m`: Extract texture components for gamma phase
  - Similar to alpha_component_extraction.m
  - Focuses on gamma phase pole figure data

- `toy_component_extraction.m`: Example script for texture component extraction
  - Demonstrates basic workflow with a toy crystal symmetry
  - Sets up for extracting multiple orientations

#### Utility Scripts
- `obj1.m` and `obj2.m`: Objective functions for computing component Orientation Distribution Functions (ODF)
  - Calculate error between observed and reconstructed ODFs
  - Support volume fraction and halfwidth parameter estimation

- `simple_orientations.m`: Reference script showing common crystallographic orientations
  - Includes predefined orientations like Cube, Brass, Copper, etc.

## Workflow

1. **Prepare Pole Figure Data**
   - Use `.xpc` files as input for pole figure data
   - Ensure correct crystal symmetry is defined

2. **Load Data**
   - Use `loadPoleFigure_beartex()` to import experimental data
   - Configure MTEX preferences for correct interpretation

3. **Extract Texture Components**
   - Run component extraction scripts (`alpha_component_extraction.m`, `gamma_component_extraction.m`)
   - Adjust parameters as needed for your specific material

4. **Generate CSV for nbragg**
   - The scripts prepare data that can be used with `nbragg.CrossSection.from_mtex()`
   - Export extracted components to a CSV format for further analysis

## MTEX Preferences
These scripts use specific MTEX preferences:
- Euler angle convention: Matthies
- X-axis direction: North
- Z-axis direction: Out of plane

## Customization
- Modify crystal symmetry in scripts to match your material
- Adjust halfwidth and volume fraction parameters in objective functions
- Add more orientation components as needed

## Troubleshooting
- Ensure MTEX toolbox is correctly installed
- Check crystal symmetry definitions
- Verify input data format

## References
- MTEX Toolbox Documentation
- Texture analysis techniques in crystallography