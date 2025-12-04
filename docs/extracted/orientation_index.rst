Crystal Orientation in nBragg
=============================

Understanding Crystal Orientation
----------------------------------

Crystal orientation is fundamental to Bragg edge neutron transmission analysis. Unlike polycrystalline materials where Bragg edges appear at fixed wavelengths regardless of orientation, single crystals produce Bragg dips whose positions depend critically on how the crystal is oriented relative to the neutron beam.

nBragg provides a straightforward framework for specifying crystal orientations while maintaining full compatibility with NCrystal's underlying physics engine.

Quick Start
-----------

The simplest way to specify orientation in nBragg:

.. code-block:: python

   import nbragg
   
   # Define crystal orientation
   # [hkl]_z: Miller indices along beam direction (z-axis)
   # [hkl]_y: Miller indices along y-axis
   
   crystal = nbragg.Crystal(
       material='Fe',
       orientation_z=[0, 0, 1],  # (001) plane perpendicular to beam
       orientation_y=[0, 1, 0]   # (010) plane along y-axis
   )

Why Orientation Matters
-----------------------

**For Single Crystals**
   The position and intensity of Bragg dips in the transmission spectrum depend on the crystal's orientation relative to the beam. Accurate orientation specification enables:
   
   - Texture analysis of materials
   - Strain/stress measurements via Bragg dip shifts
   - Crystallographic phase identification
   - Quantitative microstructure characterization

**For Polycrystals**
   While Bragg edges are orientation-independent, understanding orientation becomes crucial when analyzing:
   
   - Textured materials with preferred orientations
   - Transition from polycrystalline to single-crystal behavior
   - Orientation distribution functions (ODFs)

Coordinate System Convention
-----------------------------

nBragg uses a **fixed laboratory coordinate system** to simplify orientation specification:

- **Incident beam**: Always along the positive z-axis (:math:`\mathbf{k}_i \parallel +\mathbf{z}`)
- **Laboratory frame**: Right-handed orthogonal system :math:`(x, y, z)`
- **Crystal frame**: Defined by two reciprocal lattice vectors

This differs from NCrystal's fully flexible approach but makes orientation specification more intuitive for neutron imaging experiments where the beam direction is typically fixed.

.. note::
   This simplification does not limit the analysis - any crystal orientation can still be described using the two-vector specification method.

Documentation Structure
-----------------------

We provide three complementary resources:

.. list-table::
   :widths: 25 50 25
   :header-rows: 1

   * - Resource
     - Best For
     - Level
   * - :doc:`crystal_orientation_guide`
     - Complete theory and detailed examples
     - Beginner to Advanced
   * - :doc:`orientation_quick_reference`
     - Quick lookup during analysis
     - Intermediate to Advanced
   * - :doc:`../examples/orientation_examples`
     - Working Python code
     - All Levels

Getting Started
---------------

**New to crystal orientations?**
   Start with the :ref:`basic-concept` section of the :doc:`crystal_orientation_guide` to understand the fundamentals.

**Experienced with crystallography?**
   Jump directly to the :doc:`orientation_quick_reference` for equations and code snippets.

**Learning by doing?**
   Download and run :doc:`../examples/orientation_examples` to see practical implementations.

Common Use Cases
----------------

1. **Standard Cubic Orientation**
   
   Most common starting point:
   
   .. code-block:: python
   
      orientation_z = [0, 0, 1]
      orientation_y = [0, 1, 0]

2. **Rotated Crystal (Texture Analysis)**
   
   For textured samples:
   
   .. code-block:: python
   
      orientation_z = [1, 1, 1]  # {111} texture
      orientation_y = [1, -1, 0]

3. **Fitting to Experimental Data**
   
   Refine orientation to match Bragg dip positions:
   
   .. code-block:: python
   
      # Start with initial guess
      # Apply small rotations φ (x-axis) and θ (y-axis)
      # Minimize difference with experimental spectrum

Key Concepts
------------

Before diving into the detailed documentation, familiarize yourself with these concepts:

**Miller Indices [hkl]**
   Describe crystal planes in reciprocal lattice space. The vector :math:`\mathbf{n}_{hkl}` is perpendicular to the (hkl) crystal plane.

**Reciprocal Space**
   Orientation is specified using reciprocal lattice vectors, not real-space crystal axes. For cubic systems, these are parallel, but for other crystal systems they differ.

**Mosaicity**
   Real crystals have imperfections - crystal planes are spread over a range of orientations described by a Gaussian distribution with FWHM :math:`\eta`.

**Bragg Circles**
   At a given wavelength, neutrons satisfying the Bragg condition form a circle in reciprocal space. The crystal orientation determines if this circle intersects the crystal planes.

Mathematical Framework
----------------------

The orientation is mathematically described by a rotation :math:`g \in SO(3)` that maps the laboratory frame to the crystal frame:

.. math::

   g(\alpha, \beta, \gamma) = R_Z(\alpha) R_Y(\beta) R_Z(\gamma)

where :math:`\alpha, \beta, \gamma` are Euler angles in Matthies convention.

In practice, you specify this rotation implicitly through the two Miller index vectors **[hkl]_z** and **[hkl]_y**, which nBragg converts internally to the appropriate rotation.

See the :doc:`crystal_orientation_guide` for complete mathematical details.

Next Steps
----------

1. **Read the theory**: :doc:`crystal_orientation_guide`
2. **Try the examples**: :doc:`../examples/orientation_examples`
3. **Keep the reference handy**: :doc:`orientation_quick_reference`
4. **Apply to your data**: See :doc:`../tutorials/bragg_dip_fitting`

Additional Resources
--------------------

- **NCrystal Documentation**: For underlying physics implementation
- **Crystallography Primers**: For Miller indices and reciprocal space basics
- **nBragg Tutorials**: For complete analysis workflows

Questions?
----------

If you encounter issues with orientation specification:

1. Check the :ref:`troubleshooting-checklist` in the quick reference
2. Verify orthogonality of your orientation vectors
3. Ensure Miller indices are appropriate for your crystal system
4. Review the examples for similar crystal systems
5. Contact the nBragg community via [your support channel]

.. toctree::
   :maxdepth: 2
   :hidden:

   crystal_orientation_guide
   orientation_quick_reference
