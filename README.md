# nbragg

[![Documentation Status](https://readthedocs.org/projects/nbragg/badge/?version=latest)](https://nbragg.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

`nbragg` is a package designed for fitting neutron Bragg edge data using NCrystal cross-sections. This tool provides a straightforward way to analyze neutron transmission through polycrystalline materials, leveraging Bragg edges to extract information on material structure and composition.

## Features

- **Flexible cross-section calculations**: Interfaces with NCrystal to fetch cross-sections for crystalline materials.
- **Built-in tools for response and background functions**: Includes predefined models for instrument response (e.g., Gaussian, exponential) and background components (polynomial functions).
- **LMFit integration**: Allows flexible, nonlinear fitting of experimental data using the powerful `lmfit` library.
- **Pythonic API**: Simple-to-use, yet flexible enough for custom modeling.
- **Plotting utilities**: Provides ready-to-use plotting functions for easy visualization of results.
- **Bragg Edge Analysis**: Perform Bragg edge fitting to extract information such as d-spacing, strain, and texture.

## Installation

From source:
```bash
git clone https://github.com/TsvikiHirsh/nbragg
cd nbragg
pip install .
```
## Usage
```python
# Import nbragg
import nbragg

# Load material from NCrystal
Fe = nbragg.CrossSection.from_material("Fe_sg229_Iron-alpha")

# Load transmission data
data = nbragg.Data.from_transmission("iron_alpha.dat")

# Define a Bragg edge model
model = nbragg.TransmissionModel(Fe, vary_response=True, vary_background=True)

# Fit data using lmfit
result = model.fit(data, emin=0.4e6, emax=1.7e6)

# Plot the fit results
result.plot()
```

- Perform Bragg edge analysis on various crystalline materials using nbragg.
- Define and customize different background and response models.
- Refer to our Jupyter notebook demo for detailed usage examples.

For more detailed examples and advanced usage, please refer to our documentation page.
## License

nbragg is licensed under the MIT License.