================
Advanced Fitting
================

This guide covers advanced fitting strategies and techniques in nbragg.

Multi-Stage Refinement
=======================

nbragg supports two refinement strategies:

Rietveld Method (Default)
--------------------------

Parameters accumulate across stages - earlier parameters remain active:

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_basic=True,
        vary_background=True,
        vary_weights=True
    )

    model.stages = {
        'basic': ['norm', 'thickness'],
        'background': 'background',
        'weights': 'weights'
    }

    # Stage 1: Fit norm, thickness
    # Stage 2: Continue with norm, thickness + background
    # Stage 3: Continue with all + weights
    result = model.fit(data, wlmin=1.0, wlmax=6.0, method='rietveld')

Staged Method
-------------

Only current stage parameters vary:

.. code-block:: python

    # Stage 1: Only norm, thickness
    # Stage 2: Only background (freeze basic)
    # Stage 3: Only weights (freeze others)
    result = model.fit(data, wlmin=1.0, wlmax=6.0, method='staged')

Custom Stage Definitions
=========================

By Parameter Names
------------------

.. code-block:: python

    model.stages = {
        'scale': ['norm', 'thickness'],
        'weights_only': ['p1', 'p2']  # Explicit parameter names
    }

By Groups
---------

.. code-block:: python

    model.stages = {
        'basic_params': 'basic',      # Predefined group
        'bg': 'background',
        'phase_fractions': 'weights'
    }

Wavelength-Specific Stages
---------------------------

.. code-block:: python

    model.stages = {
        'low_wl': 'basic wlmin=1.0 wlmax=3.0',
        'high_wl': 'weights wlmin=3.0 wlmax=6.0'
    }

See Also
========

- :doc:`basic_usage` - Fundamentals
- :doc:`model_parameters` - Parameter control
- :doc:`/api/models` - API reference
