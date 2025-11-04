# nbragg
<p align="center">
<img src="https://raw.githubusercontent.com/TsvikiHirsh/nbragg/refs/heads/master/docs/source/_static/nbragg_logo.png" alt="nbragg logo" width="200"/>
</p>

[![Documentation Status](https://readthedocs.org/projects/nbragg/badge/?version=latest)](https://nbragg.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-version]][pypi-link]

nbragg is a package designed for fitting neutron Bragg edge data using [NCrystal](https://github.com/mctools/ncrystal) cross-sections. This tool provides a straightforward way to analyze neutron transmission through polycrystalline materials, leveraging Bragg edges to extract information on material structure and composition.

## Features

- **Flexible Cross-Section Calculations**: Interfaces with NCrystal to fetch cross-sections for crystalline materials.
- **SANS Modeling**: Built-in support for Small Angle Neutron Scattering (SANS) using hard-sphere models for samples with nanoscale features.
- **Extinction Effects**: Support for primary and secondary extinction modeling for large crystallites and thick samples.
- **Built-In Tools for Response and Background Functions**: Includes predefined models for instrument response (e.g., Jorgensen, square) and background components (polynomial functions).
- **LMFit Integration**: Allows flexible, nonlinear fitting of experimental data using the powerful lmfit library.
- **Rietveld-type analysis**: Enables iterative, parametric refinement of Bragg edge data using the Rietveld method, accumulating parameters across stages for robust fitting.
- **Pythonic API**: Simple-to-use, yet flexible enough for custom modeling.
- **Plotting Utilities**: Provides ready-to-use plotting functions for easy visualization of results.
- **Bragg Edge Analysis**: Perform Bragg edge fitting to extract information such as d-spacing, strain, and texture.

## Installation

### Basic Installation

To install the base package:

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

## Usage

Here's a quick example to get started:

```python
import nbragg

data = nbragg.Data.from_transmission("iron_powder.csv") # read data
xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat") # define sample
model = nbragg.TransmissionModel(xs, vary_background=True, vary_response=True) # define model
result = model.fit(data) # perform fit
result.plot() # plot results
```

![Fit Results](https://raw.githubusercontent.com/TsvikiHirsh/nbragg/refs/heads/master/notebooks/fit_results.png)

## Tutorials and Documentation

For more detailed examples and advanced usage, including custom stage definitions and Rietveld fitting, please refer to our [documentation page](https://nbragg.readthedocs.io) and check out the updated [Jupyter notebook tutorial](https://github.com/TsvikiHirsh/nbragg/blob/master/notebooks/nbragg_tutorial.ipynb).

## License

nbragg is licensed under the MIT License.

[pypi-version]: https://img.shields.io/pypi/v/nbragg.svg
[pypi-link]: https://pypi.org/project/nbragg/
[pypi-platforms]: https://img.shields.io/badge/platforms-linux%20%7C%20osx%20%7C%20windows-blue.svg