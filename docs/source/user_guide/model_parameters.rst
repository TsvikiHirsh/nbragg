================
Model Parameters
================

This guide covers all parameter control options in ``TransmissionModel``.

Overview
========

The ``TransmissionModel`` class provides fine-grained control over which parameters
vary during fitting through ``vary_*`` arguments. This allows you to:

- Fix certain parameters while varying others
- Build multi-stage fits with precise control
- Prevent overfitting by limiting free parameters
- Match your physical constraints

Parameter Categories
====================

Basic Parameters (``vary_basic``)
----------------------------------

Controls **thickness** and **norm** parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_basic=True)

**Parameters controlled:**
    - ``thickness``: Sample thickness (cm)
    - ``norm``: Normalization factor

**Default:** ``None`` (basic stage included if other stages exist)

**Use ``vary_basic=False`` when:**
    - You know the exact sample thickness
    - You want to fix the normalization
    - Fitting other parameters first

**Example:**

.. code-block:: python

    # Fix thickness and norm, fit only weights
    model = nbragg.TransmissionModel(
        xs,
        vary_basic=False,
        vary_weights=True
    )

    # Set the fixed values
    model.params['thickness'].set(value=2.0, min=0.1, max=10.0)
    model.params['norm'].set(vary=False)

    result = model.fit(data)

.. note::
    The ``temp`` (temperature) parameter is always fixed by default (``vary=False``).

Phase Weights (``vary_weights``)
---------------------------------

Controls phase fraction parameters for multi-phase materials.

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_weights=True)

**Parameters controlled:**
    - ``p1, p2, ..., pN`` for N phases (logit-transformed)
    - Derived phase fractions (e.g., ``alpha``, ``gamma``)

**How it works:**

For N phases, nbragg creates N-1 free parameters (``p1`` through ``p_{N-1}``)
and derives all phase fractions to ensure they sum to 1:

.. math::

    w_i = \\frac{e^{p_i}}{1 + \\sum_{j=1}^{N-1} e^{p_j}}

    w_N = \\frac{1}{1 + \\sum_{j=1}^{N-1} e^{p_j}}

**Example:**

.. code-block:: python

    # Two-phase material
    xs = nbragg.CrossSection(alpha="Fe_sg229_Iron-alpha.ncmat") * 0.7 + \
         nbragg.CrossSection(gamma="Fe_sg225_Iron-gamma.ncmat") * 0.3

    model = nbragg.TransmissionModel(xs, vary_weights=True)

    # Initial weights: alpha=0.7, gamma=0.3
    # Fit will optimize these fractions
    result = model.fit(data)

    # Check final fractions
    print(f"Alpha: {result.params['alpha'].value:.3f}")
    print(f"Gamma: {result.params['gamma'].value:.3f}")

Background (``vary_background``)
---------------------------------

Controls polynomial background parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_background=True,
        background="polynomial3"  # 0th, 1st, 2nd order terms
    )

**Parameters controlled:**
    - ``bg0``, ``bg1``, ``bg2``, ... (depending on polynomial order)

**Background types:**
    - ``"polynomial0"``: Constant offset
    - ``"polynomial1"``: Linear
    - ``"polynomial2"``: Quadratic
    - ``"polynomial3"``: Cubic (default)

Instrument Response (``vary_response``)
----------------------------------------

Controls instrument resolution function parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_response=True,
        response="jorgensen"
    )

**Parameters controlled (Jorgensen function):**
    - ``α0``: Shape parameter alpha
    - ``β0``: Shape parameter beta

**Response types:**
    - ``"jorgensen"``: Standard Jorgensen function (default)
    - ``"gaussian"``: Simple Gaussian

Orientation (``vary_orientation``)
-----------------------------------

Controls crystal orientation parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_orientation=True)

**Parameters controlled (per phase):**
    - ``θ_phase``: Theta rotation (degrees)
    - ``ϕ_phase``: Phi rotation (degrees)
    - ``η_phase``: Mosaicity (degrees)

**Initial values:**
    Parameters are automatically initialized from the material dictionary:

    - ``θ`` from ``material['theta']``
    - ``ϕ`` from ``material['phi']``
    - ``η`` from ``material['mos']`` (mosaicity)

**Example:**

.. code-block:: python

    # Oriented material with known mosaicity
    materials = {
        'phase1': {
            'mat': 'Fe_sg229_Iron-alpha.ncmat',
            'mos': 22.5,  # Initial mosaicity
            'theta': 10.0,
            'phi': 15.0,
            'dir1': [1, 0, 0],
            'dir2': [0, 1, 0],
            'weight': 1.0
        }
    }

    xs = nbragg.CrossSection(materials)

    # Model will initialize η_phase1 = 22.5
    model = nbragg.TransmissionModel(xs, vary_orientation=False)

    # Check initial values
    print(f"η: {model.params['η_phase1'].value}")  # 22.5
    print(f"θ: {model.params['θ_phase1'].value}")  # 10.0

See :doc:`orientation_index` for complete orientation documentation.

TOF (``vary_tof``)
------------------

Controls time-of-flight correction parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_tof=True,
        tof_length=9.0  # meters
    )

**Parameters controlled:**
    - ``L0``: Effective flight path length
    - ``t0``: Time offset

Lattice (``vary_lattice``)
---------------------------

Controls lattice parameters for phases.

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_lattice=True)

**Parameters controlled:**
    - ``a``, ``b``, ``c``: Lattice constants
    - ``a_phase``, ``b_phase``, ``c_phase``: Per-phase lattice

Extinction (``vary_extinction``)
---------------------------------

Controls extinction parameters (requires CrysExtn plugin).

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_extinction=True)

**Parameters controlled:**
    - ``ext_phase``: Extinction parameter per phase

SANS (``vary_sans``)
--------------------

Controls SANS hard-sphere radius parameters.

.. code-block:: python

    model = nbragg.TransmissionModel(xs, vary_sans=True)

**Parameters controlled:**
    - ``sans`` or ``sans_phase``: Hard-sphere radius (Angstroms)

.. note::
    Requires ``spglib`` package: ``pip install spglib``

Common Patterns
===============

Fix Everything Except Weights
------------------------------

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_basic=False,
        vary_weights=True,
        vary_background=False,
        vary_response=False
    )

    # Set known values
    model.params['thickness'].set(value=2.0)
    model.params['norm'].set(value=1.0, vary=False)

Sequential Refinement
---------------------

.. code-block:: python

    # Stage 1: Fit basic parameters
    model1 = nbragg.TransmissionModel(xs, vary_basic=True)
    result1 = model1.fit(data)

    # Stage 2: Fix basic, fit weights
    model2 = nbragg.TransmissionModel(xs, vary_basic=False, vary_weights=True)
    model2.params.update(result1.params)  # Use previous values
    result2 = model2.fit(data)

Full Rietveld Refinement
-------------------------

.. code-block:: python

    model = nbragg.TransmissionModel(
        xs,
        vary_basic=True,
        vary_background=True,
        vary_weights=True,
        vary_response=True
    )

    # Define stage order
    model.stages = {
        'basic': ['norm', 'thickness'],
        'background': 'background',
        'weights': 'weights',
        'response': 'response'
    }

    result = model.fit(data, wlmin=1.0, wlmax=6.0)

Parameter Bounds and Constraints
=================================

Setting Bounds
--------------

.. code-block:: python

    # After model creation
    model.params['thickness'].set(min=0.1, max=10.0)
    model.params['norm'].set(min=0.5, max=2.0)

Fixing Parameters
-----------------

.. code-block:: python

    # Fix at specific value
    model.params['thickness'].set(value=2.0, vary=False)

    # Fix temperature (always recommended)
    model.params['temp'].set(vary=False)

Parameter Expressions
---------------------

Some parameters are derived from others using expressions:

.. code-block:: python

    # For N=2 phases:
    # alpha = exp(p1) / (1 + exp(p1))
    # gamma = 1 / (1 + exp(p1))

    # These are automatically created and maintained

Troubleshooting
===============

**Parameters not varying?**
    Check that ``vary_*`` flag is ``True`` and parameter has ``vary=True``

**Fit diverging?**
    - Tighten parameter bounds
    - Fix some parameters initially
    - Use staged refinement

**Wrong initial values?**
    - For orientation parameters, check material dictionary
    - Set explicitly: ``model.params['param'].set(value=...)``

**Temperature varying unexpectedly?**
    Temperature is fixed by default. If it varies, something changed it.

See Also
========

- :doc:`basic_usage` - Getting started
- :doc:`advanced_fitting` - Complex fitting strategies
- :doc:`orientation_index` - Orientation parameter details
- :doc:`/api/models` - Complete API reference
