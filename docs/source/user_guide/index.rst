===========
User Guide
===========

This comprehensive guide covers all aspects of using nbragg for neutron Bragg edge transmission analysis.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   basic_usage
   model_parameters
   orientation_index
   crystal_orientation_guide
   orientation_quick_reference
   mtex_integration
   advanced_fitting

Quick Links
-----------

**New to nbragg?**
   Start with :doc:`basic_usage` to learn the fundamentals.

**Working with oriented materials?**
   See the :doc:`orientation_index` for a complete guide to crystal orientations.

**Importing MTEX data?**
   Check out :doc:`mtex_integration` for seamless MTEX workflow integration.

**Need advanced fitting control?**
   Explore :doc:`model_parameters` and :doc:`advanced_fitting`.

Overview
--------

Core Concepts
~~~~~~~~~~~~~

**Cross Sections**
   nbragg uses NCrystal to calculate neutron cross sections for crystalline materials.
   Learn how to define materials, handle multi-phase samples, and work with oriented crystals.

**Transmission Models**
   The ``TransmissionModel`` class provides flexible fitting capabilities with multi-stage
   refinement strategies including true Rietveld refinement.

**Parameter Control**
   Fine-grained control over which parameters vary during fitting, including basic parameters
   (thickness, norm), weights, orientations, and more.

Key Features
~~~~~~~~~~~~

- **Multi-phase Analysis**: Model powder and textured samples with multiple phases
- **Oriented Materials**: Full support for single crystal and textured polycrystals
- **MTEX Integration**: Import orientation distributions from MTEX for texture analysis
- **Flexible Fitting**: Rietveld and staged refinement strategies
- **Parameter Management**: Control all aspects of your fit with ``vary_*`` parameters

What's New
----------

Recent additions to nbragg include:

**Model Parameter Control**
   - New ``vary_basic`` parameter to control thickness and norm fitting
   - Automatic initialization of orientation parameters from material properties
   - Temperature parameter now fixed by default

**MTEX Integration Improvements**
   - ``from_mtex()`` now accepts material as string or dictionary
   - Automatic ``short_name`` generation from material filenames
   - Better handling of powder phases in oriented materials

**Enhanced Documentation**
   - Comprehensive crystal orientation guide
   - Quick reference cards for common operations
   - Executable examples for all major features

Getting Help
------------

- :doc:`/quickstart` - Get started quickly
- :doc:`/api/index` - Complete API reference
- :doc:`/examples/iron_powder` - Working examples
- `GitHub Issues <https://github.com/TsvikiHirsh/nbragg/issues>`_ - Report bugs or request features
