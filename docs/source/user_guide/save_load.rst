=============================
Saving and Loading Fits
=============================

nbragg provides a flexible and robust system for saving and loading both models and fit results. This allows you to:

- Save fit results for later analysis
- Share fitted models with collaborators
- Resume analysis sessions
- Archive successful fits
- Reconstruct models with fitted parameters

The save/load system avoids pickle-related issues (especially with NCrystal's ctypes) by using JSON serialization, making files human-readable and portable.

Quick Start
===========

Basic Workflow
--------------

.. code-block:: python

   from nbragg import CrossSection, TransmissionModel, materials

   # Create and fit a model
   xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
   model = TransmissionModel(xs, vary_background=True)
   result = model.fit(data)

   # Save the fit result
   result.save("my_fit.json")

   # Later, load it back
   loaded_model = TransmissionModel.load("my_fit.json")
   loaded_model.result.plot()  # Access the result with all methods

Key Features
============

Result Objects Have ``save()`` Method
--------------------------------------

All fit results now have a ``save()`` method that makes saving straightforward:

.. code-block:: python

   result = model.fit(data)
   result.save("my_result.json")

This automatically creates two files:

- ``my_result.json`` - Contains fit statistics and parameters
- ``my_result_model.json`` - Contains the model configuration

Smart File Loading
------------------

``TransmissionModel.load()`` automatically detects whether you're loading a model or a result file:

.. code-block:: python

   # Load a model configuration
   model = TransmissionModel.load("my_model.json")

   # Load a fit result (automatically loads the model too)
   model = TransmissionModel.load("my_result.json")
   model.result.redchi  # Access fit statistics

Initialize from File
--------------------

You can initialize a ``TransmissionModel`` directly from a saved file:

.. code-block:: python

   # From a model file
   model = TransmissionModel("my_model.json")

   # From a result file
   model = TransmissionModel("my_result.json")
   model.result.plot()  # Has the result attached

Complete Examples
=================

Saving and Loading Models
--------------------------

Save just the model configuration (without fit results):

.. code-block:: python

   from nbragg import CrossSection, TransmissionModel, materials

   # Create a model with specific settings
   xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
   model = TransmissionModel(
       xs,
       vary_background=True,
       vary_response=True,
       tof_length=12.5
   )

   # Define custom stages
   model.stages = {
       'background': 'background',
       'scale': ['norm', 'thickness'],
       'response': 'response'
   }

   # Save the model
   model.save("my_model.json")

   # Load it back
   loaded_model = TransmissionModel.load("my_model.json")
   # or
   loaded_model = TransmissionModel("my_model.json")

   # Stages are preserved
   print(loaded_model.stages)
   # TOF length is preserved
   print(loaded_model.tof_length)  # 12.5

Saving and Loading Fit Results
-------------------------------

Save complete fit results with all statistics and parameters:

.. code-block:: python

   # Fit the model
   result = model.fit(data, wlmin=1.0, wlmax=5.0)

   # Save using result.save()
   result.save("my_fit.json")

   # Load the result
   loaded_model = TransmissionModel.load("my_fit.json")

   # Access all fit information
   print(f"Reduced chi-square: {loaded_model.result.redchi}")
   print(f"AIC: {loaded_model.result.aic}")
   print(f"BIC: {loaded_model.result.bic}")

   # Access fitted parameters
   for param_name, param in loaded_model.result.params.items():
       if param.vary:
           print(f"{param_name}: {param.value:.6f} ± {param.stderr:.6f}")

   # Use all result methods
   loaded_model.result.plot()
   loaded_model.result.plot_total_xs()

Re-saving Loaded Results
-------------------------

Loaded results can be saved again, useful for archiving or sharing:

.. code-block:: python

   # Load a result
   model = TransmissionModel.load("original_fit.json")

   # Save it with a new name
   model.result.save("archived_fit.json")

   # Or modify parameters and save
   model.params['thickness'].value = 2.5
   new_result = model.fit(data)
   new_result.save("modified_fit.json")

Multi-Phase Materials
---------------------

Save/load works seamlessly with multi-phase materials:

.. code-block:: python

   # Create a multi-phase cross section
   xs = CrossSection({
       'alpha': {
           'mat': 'Fe_sg229_Iron-alpha.ncmat',
           'weight': 0.7
       },
       'gamma': {
           'mat': 'Fe_sg225_Iron-gamma.ncmat',
           'weight': 0.3
       }
   })

   model = TransmissionModel(xs, vary_background=True, vary_weights=True)
   result = model.fit(data)

   # Save the result
   result.save("multiphase_fit.json")

   # Load it back - all phases are preserved
   loaded_model = TransmissionModel.load("multiphase_fit.json")
   print(loaded_model._materials.keys())  # Both phases present

Advanced Usage
==============

Workflow Integration
--------------------

Integrate save/load into your analysis workflow:

.. code-block:: python

   def fit_and_save(data, sample_name):
       """Fit data and save result with organized naming."""
       xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
       model = TransmissionModel(xs, vary_background=True)

       result = model.fit(data)

       # Save with descriptive name
       filename = f"fits/{sample_name}_{result.redchi:.4f}.json"
       result.save(filename)

       return result

   def load_best_fit(sample_name):
       """Load the best fit for a sample."""
       import glob
       import json

       # Find all fits for this sample
       pattern = f"fits/{sample_name}_*.json"
       fit_files = glob.glob(pattern)

       # Load and compare
       best_redchi = float('inf')
       best_model = None

       for fit_file in fit_files:
           model = TransmissionModel.load(fit_file)
           if model.result.redchi < best_redchi:
               best_redchi = model.result.redchi
               best_model = model

       return best_model

Comparing Fits
--------------

Load multiple fits and compare them:

.. code-block:: python

   # Load several fits
   models = [
       TransmissionModel.load("fit1.json"),
       TransmissionModel.load("fit2.json"),
       TransmissionModel.load("fit3.json")
   ]

   # Compare statistics
   for i, model in enumerate(models, 1):
       print(f"Fit {i}:")
       print(f"  Reduced χ²: {model.result.redchi:.4f}")
       print(f"  AIC: {model.result.aic:.2f}")
       print(f"  BIC: {model.result.bic:.2f}")

   # Find the best fit by AIC
   best_model = min(models, key=lambda m: m.result.aic)
   print(f"\nBest model has AIC = {best_model.result.aic:.2f}")

Parameter Evolution Tracking
----------------------------

Track how parameters change across multiple fits:

.. code-block:: python

   import pandas as pd

   fit_files = ["fit_100K.json", "fit_200K.json", "fit_300K.json"]
   temperatures = [100, 200, 300]

   # Collect fitted parameters
   data = []
   for temp, fit_file in zip(temperatures, fit_files):
       model = TransmissionModel.load(fit_file)
       data.append({
           'temperature': temp,
           'thickness': model.result.params['thickness'].value,
           'norm': model.result.params['norm'].value,
           'redchi': model.result.redchi
       })

   df = pd.DataFrame(data)
   print(df)

What Gets Saved?
================

Model Files
-----------

When you save a model, the following information is stored:

- **Material specifications**: All phase definitions and properties
- **Parameters**: All parameter values, bounds, vary flags, and expressions
- **Configuration**: TOF length, response type, background type
- **Stages**: Fitting stage definitions
- **Cross-section metadata**: Name, total weight, extinction parameters

Result Files
------------

When you save a result, the following information is stored:

- **Fitted parameters**: Final parameter values with uncertainties
- **Initial parameters**: Starting values before fitting
- **Fit statistics**: χ², reduced χ², AIC, BIC
- **Fit details**: Number of function evaluations, data points, free parameters
- **Success status**: Whether the fit converged
- **Model configuration**: Complete model state (saved in separate _model.json file)

File Format
-----------

Files are saved as human-readable JSON:

.. code-block:: json

   {
     "version": "1.0",
     "class": "ModelResult",
     "params": "...",
     "init_params": "...",
     "success": true,
     "chisqr": 12.345,
     "redchi": 0.678,
     "aic": -123.45,
     "bic": -98.76,
     "nfev": 42,
     "nvarys": 5,
     "ndata": 100
   }

Best Practices
==============

Naming Conventions
------------------

Use descriptive filenames that include:

- Sample identifier
- Temperature or other conditions
- Date/time stamp (optional)
- Quality metric (optional)

.. code-block:: python

   # Good naming examples
   result.save(f"iron_powder_RT_{result.redchi:.4f}.json")
   result.save(f"steel_100K_2024-01-15.json")
   result.save(f"{sample_id}_fit.json")

Organization
------------

Organize your fits in a directory structure:

.. code-block:: text

   project/
   ├── models/
   │   ├── iron_powder_model.json
   │   └── steel_texture_model.json
   ├── fits/
   │   ├── iron_powder_RT_fit.json
   │   ├── iron_powder_100K_fit.json
   │   └── steel_texture_fit.json
   └── analysis/
       └── compare_fits.py

Version Control
---------------

JSON files work well with version control:

.. code-block:: bash

   # Add fits to git
   git add fits/*.json
   git commit -m "Add temperature series fits"

Archiving
---------

Archive important fits with metadata:

.. code-block:: python

   import json
   from datetime import datetime

   # Save the fit
   result.save("fit.json")

   # Create metadata file
   metadata = {
       'date': datetime.now().isoformat(),
       'sample': 'Iron powder RT',
       'operator': 'John Doe',
       'notes': 'High quality fit, use for publication',
       'redchi': result.redchi,
       'fit_file': 'fit.json'
   }

   with open("fit_metadata.json", 'w') as f:
       json.dump(metadata, f, indent=2)

Troubleshooting
===============

Missing Model File
------------------

If you try to load a result file without its corresponding model file:

.. code-block:: python

   # This will fail if my_result_model.json is missing
   model = TransmissionModel.load("my_result.json")
   # FileNotFoundError: Model file my_result_model.json not found

**Solution**: Ensure both files are in the same directory.

Incompatible Versions
---------------------

If loading an older version:

.. code-block:: python

   # Warning: Loading file saved with version 0.9, current version is 1.0
   model = TransmissionModel.load("old_fit.json")

**Solution**: Re-fit and save with the current version if issues arise.

Large File Sizes
----------------

For materials with many phases or parameters, files can be large.

**Solution**: Use compression:

.. code-block:: python

   import gzip
   import json

   # Save compressed
   with gzip.open("fit.json.gz", 'wt') as f:
       json.dump(state, f)

   # Load compressed
   with gzip.open("fit.json.gz", 'rt') as f:
       state = json.load(f)

API Reference
=============

``result.save(filename)``
-------------------------

Save a fit result to a JSON file.

**Parameters:**
   - ``filename`` (str): Path to the output JSON file

**Creates:**
   - ``filename``: Fit result data
   - ``filename.replace('.json', '_model.json')``: Model configuration

**Example:**

.. code-block:: python

   result.save("my_fit.json")

``TransmissionModel.save(filename)``
------------------------------------

Save the model configuration to a JSON file.

**Parameters:**
   - ``filename`` (str): Path to the output JSON file

**Example:**

.. code-block:: python

   model.save("my_model.json")

``TransmissionModel.load(filename)``
------------------------------------

Load a model or result from a JSON file.

**Parameters:**
   - ``filename`` (str): Path to the input JSON file (model or result)

**Returns:**
   - ``TransmissionModel``: The reconstructed model
   - If loading a result, ``model.result`` contains the loaded fit result

**Example:**

.. code-block:: python

   # Load a model
   model = TransmissionModel.load("my_model.json")

   # Load a result
   model = TransmissionModel.load("my_result.json")
   print(model.result.redchi)

``TransmissionModel(filename)``
-------------------------------

Initialize a TransmissionModel from a saved file.

**Parameters:**
   - First argument can be either a CrossSection object or a filename string

**Example:**

.. code-block:: python

   # Normal initialization
   model = TransmissionModel(cross_section, vary_background=True)

   # Initialize from file
   model = TransmissionModel("my_model.json")
   model = TransmissionModel("my_result.json")

See Also
========

- :doc:`basic_usage` - Getting started with nbragg
- :doc:`model_parameters` - Understanding model parameters
- :doc:`advanced_fitting` - Advanced fitting strategies
- :doc:`/api/index` - Complete API reference
