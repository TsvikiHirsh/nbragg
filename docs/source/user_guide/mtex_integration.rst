================
MTEX Integration
================

This guide shows how to seamlessly integrate MTEX orientation data into nbragg for texture analysis.

Overview
========

MTEX is a powerful MATLAB toolbox for analyzing crystallographic textures. nbragg provides
direct import of MTEX orientation distributions for Bragg edge texture analysis.

**Key features:**

- Import orientation distributions from MTEX CSV files
- Automatic handling of Euler angles and crystal directions
- Built-in support for powder phases
- Flexible material specification
- Automatic phase naming

Quick Start
===========

Basic Import
------------

.. code-block:: python

    import nbragg

    # Import MTEX data - all arguments optional!
    xs = nbragg.CrossSection.from_mtex(
        "orientation_data.csv",
        material="Fe_sg225_Iron-gamma.ncmat"
    )

    # That's it! Ready to model
    model = nbragg.TransmissionModel(xs, vary_weights=True)
    result = model.fit(data)

Material Specification
======================

nbragg accepts materials in three flexible ways:

Method 1: Material Dictionary (Original)
-----------------------------------------

.. code-block:: python

    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material=nbragg.materials["Fe_sg225_Iron-gamma.ncmat"],
        short_name="gamma"
    )

Method 2: Material String (New!)
---------------------------------

.. code-block:: python

    # Just pass the filename as a string
    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        short_name="gamma"
    )

Method 3: Auto-Generated Names (New!)
--------------------------------------

.. code-block:: python

    # Omit short_name - automatically extracts from material
    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat"
    )
    # Creates phases: Iron-gamma0, Iron-gamma1, ..., Iron-gamma_powder

CSV Format
==========

Required Columns
----------------

Your MTEX CSV must contain these columns:

.. code-block:: text

    alpha_mtex, beta_mtex, gamma_mtex  # Euler angles (degrees)
    volume_mtex                        # Volume fraction for each orientation
    fwhm                              # Full width at half maximum (mosaicity, degrees)
    xh, xk, xl                        # First crystal direction components
    yh, yk, yl                        # Second crystal direction components

Example CSV
-----------

.. code-block:: text

    alpha_mtex,beta_mtex,gamma_mtex,volume_mtex,xh,xk,xl,yh,yk,yl,fwhm
    289.73,35.96,93.54,0.0897,2.643,-0.948,0.568,1.100,2.119,-1.584,12.59
    347.56,45.15,357.89,0.0873,1.949,0.689,1.983,-0.538,2.779,-0.438,17.83
    ...

Exporting from MTEX
--------------------

Use these MATLAB commands in MTEX:

.. code-block:: matlab

    % After MTEX orientation analysis
    odf = ... % your ODF
    % Calculate orientation components
    [ori, vol] = discreteSample(odf, N);

    % Export to CSV
    T = table(ori.alpha*180/pi, ori.beta*180/pi, ori.gamma*180/pi, vol, ...
              'VariableNames', {'alpha_mtex','beta_mtex','gamma_mtex','volume_mtex'});
    writetable(T, 'orientations.csv');

See the nbragg repository for complete MTEX export scripts.

Orientation Representation
==========================

Automatic Conversion
--------------------

nbragg automatically converts MTEX orientations to NCrystal format:

1. **Euler angles** → Rotation matrices
2. **Crystal directions** → Normalized ``dir1`` and ``dir2``
3. **FWHM** → Mosaicity (``mos`` parameter)

The conversion handles:

- Coordinate system differences between MTEX and NCrystal
- Vector normalization
- Volume fraction normalization

Powder Phase
------------

By default, ``from_mtex`` adds a powder phase for the remaining volume:

.. code-block:: python

    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        powder_phase=True  # Default
    )

If your orientation volumes sum to 0.7, a powder phase with weight 0.3 is added automatically.

To disable:

.. code-block:: python

    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        powder_phase=False
    )

Phase Naming
============

Auto-Generated Names
--------------------

When ``short_name`` is not provided, nbragg extracts it from the material:

.. code-block:: python

    # Material: "Fe_sg225_Iron-gamma.ncmat"
    # Extracts: "Iron-gamma"
    # Creates: Iron-gamma0, Iron-gamma1, ..., Iron-gamma_powder

Custom Short Names
------------------

Use custom prefixes for clarity:

.. code-block:: python

    xs = nbragg.CrossSection.from_mtex(
        "alpha_orientations.csv",
        material="Fe_sg229_Iron-alpha.ncmat",
        short_name="α"
    )
    # Creates: α0, α1, α2, ..., α_powder

    xs = nbragg.CrossSection.from_mtex(
        "gamma_orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        short_name="γ"
    )
    # Creates: γ0, γ1, γ2, ..., γ_powder

Complete Example
================

Steel Texture Analysis
----------------------

.. code-block:: python

    import nbragg

    # 1. Load transmission data
    data = nbragg.Data.from_counts(
        "steel_sample.csv",
        "openbeam.csv"
    )

    # 2. Import alpha-iron texture from MTEX
    xs_alpha = nbragg.CrossSection.from_mtex(
        "alpha_texture.csv",
        material="Fe_sg229_Iron-alpha.ncmat",
        short_name="α"
    )

    # 3. Import gamma-iron texture from MTEX
    xs_gamma = nbragg.CrossSection.from_mtex(
        "gamma_texture.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        short_name="γ"
    )

    # 4. Combine phases
    xs_total = xs_alpha * 0.6 + xs_gamma * 0.4

    # 5. Create model
    model = nbragg.TransmissionModel(
        xs_total,
        vary_weights=True,      # Fit phase fractions
        vary_orientation=False,  # Keep MTEX orientations fixed
        vary_background=True,
        vary_response=True
    )

    # 6. Check initial mosaicity values from MTEX
    for phase in ['α0', 'α1', 'γ0', 'γ1']:
        if f'η_{phase}' in model.params:
            print(f"{phase} mosaicity: {model.params[f'η_{phase}'].value:.1f}°")

    # 7. Fit
    result = model.fit(data, wlmin=1.0, wlmax=6.0)

    # 8. Analyze results
    print(result.fit_report())
    result.plot()

Advanced Usage
==============

Custom Material Properties
--------------------------

Pass custom properties in the material dict:

.. code-block:: python

    custom_material = {
        'mat': 'Fe_sg225_Iron-gamma.ncmat',
        'temp': 350.0,  # Custom temperature
        'weight': 1.0
    }

    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material=custom_material,
        short_name="gamma_350K"
    )

Multiple Texture Components
----------------------------

Combine different texture components:

.. code-block:: python

    # Sharp texture (small FWHM)
    xs_sharp = nbragg.CrossSection.from_mtex(
        "sharp_texture.csv",
        material="Fe_sg229_Iron-alpha.ncmat",
        short_name="α_sharp"
    )

    # Weak texture (large FWHM)
    xs_weak = nbragg.CrossSection.from_mtex(
        "weak_texture.csv",
        material="Fe_sg229_Iron-alpha.ncmat",
        short_name="α_weak"
    )

    # Combine with weights
    xs = xs_sharp * 0.3 + xs_weak * 0.7

Accessing Orientation Data
---------------------------

After import, orientation data is in the materials dictionary:

.. code-block:: python

    xs = nbragg.CrossSection.from_mtex(
        "orientations.csv",
        material="Fe_sg225_Iron-gamma.ncmat",
        short_name="gamma"
    )

    # Access first oriented phase
    phase0 = xs.materials['gamma0']
    print(f"Mosaicity: {phase0['mos']}")
    print(f"Direction 1: {phase0['dir1']}")
    print(f"Direction 2: {phase0['dir2']}")
    print(f"Weight: {phase0['weight']}")

    # Check powder phase
    powder = xs.materials['gamma_powder']
    print(f"Powder weight: {powder['weight']}")

Troubleshooting
===============

Missing Columns
---------------

**Error:** ``KeyError: Could not find column``

**Solution:** Ensure CSV has all required columns. Check column names match exactly.

Coordinate System Issues
------------------------

**Problem:** Orientations don't match expected directions

**Solution:** Verify MTEX export uses the same coordinate conventions. The nbragg
MTEX scripts handle this automatically.

Volume Normalization
--------------------

**Problem:** Volumes don't sum to 1

**Solution:** nbragg automatically normalizes volumes. If they sum > 1, they're scaled down.

Weight Balancing
----------------

**Problem:** Powder phase weight is negative or very large

**Solution:** Check that MTEX volumes are reasonable. Typically should sum to < 1.

Integration Checklist
======================

Before using MTEX data in production:

1. ✅ Verify CSV format matches requirements
2. ✅ Test import: ``xs = CrossSection.from_mtex(...)``
3. ✅ Check phase weights sum to 1
4. ✅ Verify orientation parameters in model
5. ✅ Compare with MTEX pole figures
6. ✅ Test fit convergence
7. ✅ Validate results physically

See Also
========

- :doc:`orientation_index` - Complete orientation guide
- :doc:`model_parameters` - Parameter control
- :doc:`advanced_fitting` - Multi-component fits
- `MTEX Website <https://mtex-toolbox.github.io/>`_ - MTEX documentation
