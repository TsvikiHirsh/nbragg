# nbragg

<p align="center">
  <img src="docs/source/_static/nbragg_logo.png" alt="nbragg logo" width="200"/>
</p>

[![Documentation Status](https://readthedocs.org/projects/nbragg/badge/?version=latest)](https://nbragg.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-version]][pypi-link]


nbragg is a package designed for fitting neutron Bragg edge data using NCrystal cross-sections. This tool provides a straightforward way to analyze neutron transmission through polycrystalline materials, leveraging Bragg edges to extract information on material structure and composition.

## Features

- **Flexible Cross-Section Calculations**: Interfaces with NCrystal to fetch cross-sections for crystalline materials.
- **Built-In Tools for Response and Background Functions**: Includes predefined models for instrument response (e.g., Jorgensen, square) and background components (polynomial functions).
- **LMFit Integration**: Allows flexible, nonlinear fitting of experimental data using the powerful lmfit library.
- **Rietveld-type analysis**: Enables iterative, parametric refinement of Bragg edge data using the Rietveld method, accumulating parameters across stages for robust fitting.
- **Pythonic API**: Simple-to-use, yet flexible enough for custom modeling.
- **Plotting Utilities**: Provides ready-to-use plotting functions for easy visualization of results.
- **Bragg Edge Analysis**: Perform Bragg edge fitting to extract information such as d-spacing, strain, and texture.

## Installation

To install use:

```bash
pip install nbragg
```

## Usage

Here's a quick example to get started:

```python
import nbragg
data = nbragg.Data.from_transmission("iron_powder.csv")                         # read data
xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")                      # define sample
model = nbragg.TransmissionModel(xs, vary_background=True, vary_response=True)  # define model
result = model.fit(data)                                                        # perform fit
result.plot()                                                                   # plot results
```

![Fit Results](notebooks/fit_results.png)

## Tutorials and Documentation

For more detailed examples and advanced usage, including custom stage definitions and Rietveld fitting, please refer to our [documentation page](https://nbragg.readthedocs.io) and check out the updated [Jupyter notebook tutorial](notebooks/nbragg_tutorial.ipynb).

## License

nbragg is licensed under the MIT License.

[pypi-version]: https://img.shields.io/pypi/v/nbragg.svg
[pypi-link]: https://pypi.org/project/nbragg/
[pypi-platforms]: https://img.shields.io/badge/platforms-linux%20%7C%20osx%20%7C%20windows-blue.svg