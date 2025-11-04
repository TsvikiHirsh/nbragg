# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2025-01-XX

### Added
- **SANS (Small Angle Neutron Scattering) support**: Added hard-sphere SANS modeling capability using NCrystal's NCMATComposer
  - New `sans` parameter in material definitions to specify hard-sphere radius (Angstroms)
  - Support for `vary_sans=True` in `TransmissionModel` to fit SANS parameters
  - Multi-phase SANS with individual parameters (`sans1`, `sans2`, etc.) for each phase
  - Automatic void phase addition (0.01 fraction) as required by NCrystal
  - Integration with Rietveld staged fitting via `'sans'` group
  - Bounds: 0-1000 Angstroms for SANS radius
- **Extinction modeling documentation**: Comprehensive documentation and examples for extinction parameters in tutorial
- Tutorial section covering SANS usage with practical examples
- `register_material()` now accepts `sans` parameter for convenience

### Fixed
- **Critical bug fix**: Virtual filename collision when combining CrossSections with different SANS values
  - Previously, combining two materials from the same base file but with different SANS parameters would result in identical cross-section values due to filename collision
  - Now creates unique virtual filenames (e.g., `iron_Fe_sg225_Iron-gamma.nbragg`, `iron_1_Fe_sg225_Iron-gamma.nbragg`)
  - Added `_original_mat` field to track original .ncmat files for NCMATComposer
  - Updated 17 tests to verify correct behavior

### Changed
- `CrossSection.materials[material]['mat']` now points to unique virtual `.nbragg` files
- `CrossSection.materials[material]['_original_mat']` stores the original `.ncmat` filename
- Simplified SANS section in tutorial notebook for better readability

### Technical Details
- Modified `_process_materials()` to preserve `_original_mat` field during material processing
- Updated `_create_virtual_materials()` to use material name prefixes for unique filenames
- Enhanced `_update_ncmat_parameters()` to correctly handle SANS with NCMATComposer
- Added `_make_sans_params()` method to `TransmissionModel`
- SANS parameters integrated into stages system with `'sans'` group support

## [0.4.0] - 2025-01-XX

### Added
- Initial release features
- CrossSection modeling
- TransmissionModel with Rietveld fitting
- Extinction parameter support
- Oriented materials
- MTEX integration
- Background and response modeling

[0.5.0]: https://github.com/TsvikiHirsh/nbragg/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/TsvikiHirsh/nbragg/releases/tag/v0.4.0
