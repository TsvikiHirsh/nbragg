===========
Basic Usage
===========

This guide covers the fundamental concepts and workflows in nbragg.

Quick Example
=============

Here's a complete example showing the basic workflow:

.. code-block:: python

    import nbragg

    # 1. Load your data
    data = nbragg.Data.from_transmission("my_transmission.csv")

    # 2. Define your material
    xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")

    # 3. Create a model
    model = nbragg.TransmissionModel(
        xs,
        vary_background=True,
        vary_response=True
    )

    # 4. Fit the data
    result = model.fit(data, wlmin=1.0, wlmax=6.0)

    # 5. Plot the results
    result.plot()

Data Loading
============

nbragg supports two main data formats:

From Transmission Data
----------------------

If you have pre-calculated transmission data (wavelength, transmission, error):

.. code-block:: python

    data = nbragg.Data.from_transmission("transmission.csv")

Expected CSV format:

.. code-block:: text

    wavelength,trans,err
    1.0,0.95,0.01
    1.1,0.94,0.01
    ...

From Counts Data
----------------

If you have raw counts (signal and open beam):

.. code-block:: python

    data = nbragg.Data.from_counts(
        signal="sample_counts.csv",
        openbeam="openbeam_counts.csv",
        L=10.0,      # Flight path length in meters
        tstep=10e-6  # Time step in seconds
    )

Expected CSV format:

.. code-block:: text

    tof,counts,err
    0,1000,31.6
    1,1050,32.4
    ...

With empty region correction for background subtraction:

.. code-block:: python

    data = nbragg.Data.from_counts(
        signal="sample_counts.csv",
        openbeam="openbeam_counts.csv",
        empty_signal="empty_sample.csv",
        empty_openbeam="empty_openbeam.csv",
        L=10.0,
        tstep=10e-6
    )

Using DataFrames Instead of Files
----------------------------------

Both ``from_counts`` and ``from_transmission`` accept pandas DataFrames in addition to file paths:

.. code-block:: python

    import pandas as pd
    import nbragg

    # Create or load a DataFrame
    df = pd.DataFrame({
        'wavelength': [1.0, 1.1, 1.2],
        'trans': [0.95, 0.94, 0.93],
        'err': [0.01, 0.01, 0.01]
    })

    data = nbragg.Data.from_transmission(df)

This is useful when data is already in memory or needs preprocessing.

Time-of-Flight Corrections
---------------------------

If you've optimized L0 and t0 parameters using ``vary_tof``, you can apply these corrections when loading data:

.. code-block:: python

    # After optimizing with vary_tof=True, you get L0 and t0 values
    data = nbragg.Data.from_counts(
        signal="sample_counts.csv",
        openbeam="openbeam_counts.csv",
        L=10.0,
        tstep=10e-6,
        L0=1.05,   # Flight path scale factor (default 1.0)
        t0=5.0     # Time offset correction (tof units)
    )

The corrections are applied as:

- **L0** is a scale factor (default 1.0):

  - L0 > 1.0: longer flight path → TOF is reduced
  - L0 < 1.0: shorter flight path → TOF is increased

- **t0** is an additive time offset in TOF units

- Formula: ``corrected_tof = tof + (1.0 - L0) * tof + t0``

Use the same L0 and t0 values you obtained from fitting with ``vary_tof=True``.

Handling NaN Values
-------------------

Remove rows with NaN values during data loading:

.. code-block:: python

    # Remove NaN values when loading
    data = nbragg.Data.from_transmission("data.csv", dropna=True)

Or remove them after loading:

.. code-block:: python

    # Load data
    data = nbragg.Data.from_transmission("data.csv")

    # Remove NaN values in place
    data.dropna(inplace=True)

    # Or create a new Data object without NaN
    clean_data = data.dropna()

Defining Materials
==================

Single Phase
------------

For a single crystalline phase:

.. code-block:: python

    xs = nbragg.CrossSection(iron="Fe_sg229_Iron-alpha.ncmat")

Using the materials registry:

.. code-block:: python

    xs = nbragg.CrossSection(
        iron=nbragg.materials["Fe_sg229_Iron-alpha.ncmat"]
    )

Multi-Phase
-----------

For samples with multiple phases:

.. code-block:: python

    xs = nbragg.CrossSection(
        alpha="Fe_sg229_Iron-alpha.ncmat",
        gamma="Fe_sg225_Iron-gamma.ncmat"
    )

With initial weight fractions:

.. code-block:: python

    xs = nbragg.CrossSection(alpha="Fe_sg229_Iron-alpha.ncmat") * 0.7 + \
         nbragg.CrossSection(gamma="Fe_sg225_Iron-gamma.ncmat") * 0.3

Creating Models
===============

The ``TransmissionModel`` class is the core of nbragg fitting.

Basic Model
-----------

Simplest model with default parameters:

.. code-block:: python

    model = nbragg.TransmissionModel(xs)

Controlling Parameters
----------------------

Use ``vary_*`` flags to control what parameters can vary during fitting:

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_basic=True,        # thickness, norm
        vary_weights=True,      # phase fractions
        vary_background=True,   # background parameters
        vary_response=True,     # instrument response
        vary_orientation=False  # keep orientations fixed
    )

See :doc:`model_parameters` for complete details.

Fitting Data
============

Basic Fit
---------

Fit over a wavelength range:

.. code-block:: python

    result = model.fit(data, wlmin=1.0, wlmax=6.0)

The fitting uses the Rietveld method by default, which accumulates parameters
across stages for optimal convergence.

Custom Stages
-------------

You can define custom fitting stages:

.. code-block:: python

    model.stages = {
        'scale': ['norm', 'thickness'],
        'background': 'background',
        'weights': 'weights',
        'response': 'response'
    }

    result = model.fit(data, wlmin=1.0, wlmax=6.0)

Or pass stages directly to fit:

.. code-block:: python

    custom_stages = {
        'basic': ['norm', 'thickness'],
        'all': 'all'
    }

    result = model.fit(data, wlmin=1.0, wlmax=6.0, stages=custom_stages)

Analyzing Results
=================

The fit result contains all the information about the fit:

Parameter Values
----------------

.. code-block:: python

    # Print fit report
    print(result.fit_report())

    # Access parameters
    thickness = result.params['thickness'].value
    thickness_err = result.params['thickness'].stderr

    # Check if fit succeeded
    if result.success:
        print(f"Chi-squared: {result.chisqr}")
        print(f"Reduced chi-squared: {result.redchi}")

Stage Summary
-------------

For multi-stage fits, view the progression:

.. code-block:: python

    # Display stage summary (HTML table)
    result.stages_summary

    # Or as plain text
    print(result.stages_summary.to_html())

Visualization
-------------

.. code-block:: python

    # Plot fit result
    result.plot()

    # Plot specific stage
    result.plot_stage_progression('thickness')

    # Plot residuals
    result.plot(plot_residuals=True)

Setting Parameter Values
=========================

You can set initial values and bounds before fitting:

.. code-block:: python

    # Set value and bounds
    model.params['thickness'].set(value=2.0, min=0.1, max=10.0)

    # Fix a parameter
    model.params['norm'].set(vary=False)

    # Set with stderr
    model.params['thickness'].set(value=2.0, vary=True)

Next Steps
==========

- Learn about :doc:`model_parameters` for fine-grained control
- Explore :doc:`orientation_index` for oriented materials
- See :doc:`mtex_integration` for texture analysis
- Check :doc:`advanced_fitting` for complex scenarios
