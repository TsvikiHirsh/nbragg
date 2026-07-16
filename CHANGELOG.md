# Changelog

All notable changes to nbragg are documented in this file.

## 1.0.0 (2026-07-16)

First stable release. 🎉

### Highlights

- nbragg is now considered stable: the public API (`Data`, `CrossSection`,
  `TransmissionModel`, `GroupedFitResult`, `Response`, `Background`,
  `materials`, `register_material`, `save_result`, `load_result`) follows
  semantic versioning from this release on.
- Requires Python 3.9+.

### Fixed

- `GroupedFitResult.summary()` now returns the documented pandas DataFrame;
  a duplicate method definition that returned/displayed HTML had shadowed it.
  The HTML view is available as the new `summary_html()` method.
- Fixed `vary_orientation=True` fits crashing with
  `NCBadInput: mos must be in range (0.0,pi/2]` for cross-sections containing
  an unoriented (powder) phase — for example those imported with
  `CrossSection.from_mtex(..., powder_phase=True)`. Orientation parameters
  are no longer applied to (or varied for) unoriented phases, and mosaicity
  parameters are now bounded to NCrystal's valid (0, 90] degree range.
- `TransmissionModel.fit()` raises a clear `ValueError` when the fit
  wavelength range (`wlmin`/`wlmax`) contains no data points, instead of an
  obscure numpy error.
- The materials cache is now written atomically, so concurrent nbragg
  processes can no longer corrupt it; a corrupt cache is silently discarded
  and regenerated instead of emitting a warning.
- `ipywidgets`/`IPython` are no longer imported at package import time —
  they are only needed for `interactive_plot()` and are now an optional
  extra: `pip install "nbragg[interactive]"`.

### Changed

- `spglib` is now a regular dependency (required by NCrystal's material
  composer for SANS material construction).
- Package description corrected to "neutron Bragg edge fitting".
- Development status classifier promoted to Production/Stable.

### Documentation

- Documentation overhauled for 1.0: all user-guide pages are reachable from
  the table of contents, broken cross-references fixed, and the API reference
  now covers the grouped-fitting and save/load modules. The docs build is
  warning-free.
- All three tutorial notebooks (getting started, Rietveld refinement,
  grouped fitting) are verified to execute end-to-end against this release,
  and a notebooks index (`notebooks/README.md`) describes the recommended
  reading order.
- Added `CITATION.cff` with citation metadata.

### CI

- The GitHub Actions workflow now runs the full test suite (previously it
  only checked that the package imports and the CrysExtn plugin loads).

## 0.8.2 (2026-06-11)

- Fixed parallel grouped fitting when the model references in-memory NCrystal
  materials (user-registered via `NC.registerInMemoryFileData` or virtual
  `.nbragg` template files): the parent process now snapshots NCrystal's
  virtual-file registry and re-registers it in each worker process.

## 0.8.1 (2026-06-10)

- `Data.from_grouped()` and `Data.from_transmission()` accept dicts/lists of
  pandas DataFrames directly.
- Added `Data.from_grouped_dataframes()` and
  `Data.from_grouped_transmission()` class methods.
- Added `jorgensen_inv` response: a wavelength-resolution-invariant Jorgensen
  response function.

## 0.8.0 (2026-05-25)

- Refactored the full Jorgensen response parameters for improved stability
  and normalization, with staged fitting support.

## 0.7.1 (2026-05-07)

- Fixed several lattice-parameter fitting bugs: `NCBadInput` errors when
  reusing a `CrossSection` across models or fitting non-cubic cells, and
  `update_params` stripping lmfit expression constraints from b/c lattice
  parameters.
- More flexible column naming in `Data.from_counts()`.

## 0.7 (2026-05-06)

- Extinction parameters use physical units, with `BC_mix` as default mode.
- Added `comp` parameter to `CrossSection` for NCrystal component selection.
- Added selective group refitting (`refit`) and query-based group loading.
- Improved grouped-result visualization (colorbar labels, axis limits,
  plotting into existing axes).
- Named orientation parameters and dynamic `epsfcn` computation for lattice
  fitting.

## 0.6 (2025-12-17)

- First release with grouped/gridded data fitting: `Data.from_grouped()`,
  parallel fitting via `ProcessPoolExecutor`, `GroupedFitResult` with
  parameter maps, save/load, and flexible group indexing.

## Earlier versions

Earlier development releases (0.1–0.5) established the core workflow:
NCrystal-based cross-sections, transmission models with lmfit, Rietveld-style
staged refinement, oriented materials, MTEX integration, SANS and extinction
support. See the [GitHub releases](https://github.com/TsvikiHirsh/nbragg/releases)
for details.
