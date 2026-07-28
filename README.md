# nbragg
<p align="center">
<img src="https://raw.githubusercontent.com/TsvikiHirsh/nbragg/refs/heads/master/docs/source/_static/nbragg_logo.png" alt="nbragg logo" width="200"/>
</p>

[![CI](https://github.com/TsvikiHirsh/nbragg/actions/workflows/test.yml/badge.svg)](https://github.com/TsvikiHirsh/nbragg/actions/workflows/test.yml)
[![Documentation Status](https://readthedocs.org/projects/nbragg/badge/?version=latest)](https://nbragg.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-version]][pypi-link]
[![Python versions](https://img.shields.io/pypi/pyversions/nbragg.svg)](https://pypi.org/project/nbragg/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

nbragg is a package designed for fitting neutron Bragg edge data using [NCrystal](https://github.com/mctools/ncrystal) cross-sections. This tool provides a straightforward way to analyze neutron transmission through polycrystalline materials, leveraging Bragg edges to extract information on material structure and composition.

## Quick Example

Fit a measured iron-powder transmission spectrum in a few lines:

```python
import nbragg

data = nbragg.Data.from_transmission("iron_powder.csv") # read data
xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat") # define sample
model = nbragg.TransmissionModel(xs, vary_background=True, vary_response=True) # define model
result = model.fit(data) # perform fit
result.plot() # plot results
```

![Fit Results](https://raw.githubusercontent.com/TsvikiHirsh/nbragg/refs/heads/master/notebooks/fit_results.png)

## Features

- **Flexible Cross-Section Calculations**: Interfaces with NCrystal to fetch cross-sections for crystalline materials.
- **Multi-Phase and Oriented Materials**: Combine phases with simple arithmetic (`xs = 0.5*alpha + 0.5*gamma`), and model textured or single-crystal samples with full orientation control.
- **Grouped/Gridded Data Fitting**: Analyze spatially-resolved or multi-sample data with support for 1D arrays, 2D grids, and named groups. Includes parallel fitting with automatic result visualization via parameter maps.
- **Rietveld-Type Analysis**: Iterative, staged refinement of Bragg edge data, accumulating parameters across stages for robust fitting.
- **MTEX Integration**: Import phase weights and orientation distributions exported from [MTEX](https://mtex-toolbox.github.io/) for texture analysis.
- **SANS Modeling**: Built-in support for Small Angle Neutron Scattering (SANS) using hard-sphere models for samples with nanoscale features.
- **Extinction Effects**: Support for primary and secondary extinction modeling for large crystallites and thick samples.
- **Built-In Response and Background Functions**: Includes predefined models for instrument response (e.g., Jorgensen, square) and background components (polynomial functions).
- **LMFit Integration**: Flexible, nonlinear fitting of experimental data using the powerful lmfit library.
- **Save and Load**: JSON-based persistence of models and fit results for reproducible analysis sessions.
- **Pythonic API**: Simple to use, yet flexible enough for custom modeling.
- **Plotting Utilities**: Ready-to-use plotting functions for easy visualization of data, cross-sections, and fit results.

## Installation

nbragg requires Python 3.9 or later.

### Basic Installation

To install the base package from PyPI:

```bash
pip install nbragg
```

### Installation with Extinction Effects

To include extinction effects in your analysis, you'll need to install the extinction plugin separately:

```bash
pip install nbragg
pip install git+https://github.com/XuShuqi7/ncplugin-CrysExtn
```

The [ncrystal-plugin-crysextn](https://github.com/XuShuqi7/ncplugin-CrysExtn) plugin provides extinction corrections for crystallographic calculations.

**Note:** The extinction plugin is only required if you plan to use extinction effects. For standard Bragg edge fitting without extinction corrections, the base installation is sufficient.

## Tutorials and Documentation

Full documentation is available at [nbragg.readthedocs.io](https://nbragg.readthedocs.io).

Three Jupyter notebook tutorials cover the main workflows:

1. [Getting started with nbragg](https://github.com/TsvikiHirsh/nbragg/blob/master/notebooks/nbragg_tutorial.ipynb) — data loading, cross-sections, model definition, fitting, oriented materials, and MTEX integration.
2. [Rietveld-type refinement](https://github.com/TsvikiHirsh/nbragg/blob/master/notebooks/Rietveld_in_nbragg_tutorial.ipynb) — staged, parametric refinement of Bragg edge data.
3. [Grouped/gridded data fitting](https://github.com/TsvikiHirsh/nbragg/blob/master/notebooks/grouped_fits_tutorial.ipynb) — spatially-resolved and multi-sample analysis with parallel fitting.

## Citing nbragg

If you use nbragg in your research, please cite the accompanying paper:

> T. Y. Hirsh, A. F. T. Leong, A. M. Long, D. D. DiJulio, S. Xu, G. Muhrer,
> T. H. Kittelmann, J. I. Marquez Damian, D. J. Savage and S. C. Vogel,
> *nbragg: A Versatile Python Tool for Bragg-Edge Transmission Analysis Using NCrystal*,
> Journal of Applied Crystallography (submitted, 2026).

BibTeX:

```bibtex
@article{nbragg2026,
  title   = {nbragg: A Versatile Python Tool for Bragg-Edge Transmission Analysis Using NCrystal},
  author  = {Hirsh, Tsviki Y. and Leong, Andrew F. T. and Long, Alexander M. and
             DiJulio, Douglas D. and Xu, Shuqi and Muhrer, G{\"u}nter and
             Kittelmann, Thomas H. and Marquez Damian, Jos{\'e} I. and
             Savage, Daniel J. and Vogel, Sven C.},
  journal = {Journal of Applied Crystallography},
  year    = {2026},
  note    = {Submitted}
}
```

Citation metadata is also provided in [CITATION.cff](CITATION.cff), and GitHub's "Cite this repository" button generates BibTeX/APA entries automatically.

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines, and use the [issue tracker](https://github.com/TsvikiHirsh/nbragg/issues) to report bugs or request features.

## License

nbragg is licensed under the [MIT License](LICENSE).

### Third-Party Dependencies

This project depends on several open-source packages with permissive licenses compatible with MIT:

- **scipy, pandas, numpy, lmfit**: BSD/BSD 3-Clause
- **setuptools, tqdm**: MIT License
- **matplotlib**: PSF License
- **ncrystal**: Apache 2.0 ([license](https://github.com/mctools/ncrystal/blob/master/LICENSE))

All dependencies allow free use, modification, and distribution.

[pypi-version]: https://img.shields.io/pypi/v/nbragg.svg
[pypi-link]: https://pypi.org/project/nbragg/
[pypi-platforms]: https://img.shields.io/badge/platforms-linux%20%7C%20osx%20%7C%20windows-blue.svg
