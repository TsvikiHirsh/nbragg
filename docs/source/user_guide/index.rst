===========
User Guide
===========

This comprehensive guide covers all aspects of using nbragg for neutron Bragg edge transmission analysis.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   basic_usage
   grouped_fitting
   model_parameters
   save_load
   orientation_index
   crystal_orientation_guide
   orientation_quick_reference
   mtex_integration
   advanced_fitting

Quick Links
-----------

**New to nbragg?**
   Start with :doc:`basic_usage` to learn the fundamentals.

**Analyzing grouped or spatially-resolved data?**
   See :doc:`grouped_fitting` for 2D grids, 1D arrays, and multi-sample analysis.

**Saving and loading fits?**
   See :doc:`save_load` to learn how to save results and resume analysis sessions.

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

- **Grouped/Gridded Data**: Analyze spatially-resolved or multi-sample datasets with parallel fitting
- **Multi-phase Analysis**: Model powder and textured samples with multiple phases
- **Oriented Materials**: Full support for single crystal and textured polycrystals
- **MTEX Integration**: Import orientation distributions from MTEX for texture analysis
- **Flexible Fitting**: Rietveld and staged refinement strategies
- **Parameter Management**: Control all aspects of your fit with ``vary_*`` parameters

What's New
----------

Recent additions to nbragg include:

**Grouped/Gridded Data Fitting**
   - ``Data.from_grouped()`` for loading spatially-resolved or multi-sample data
   - Parallel fitting with ``n_jobs`` parameter
   - ``plot_parameter_map()`` with auto-detection of plot type (heatmap/line/bar)
   - Flexible indexing: tuples ``(0,0)``, strings ``"(0,0)"``, integers, or names
   - ``GroupedFitResult`` with ``save()``, ``load()``, and ``fit_report()`` methods
   - See :doc:`grouped_fitting` for complete documentation

**Save and Load Functionality**
   - ``result.save()`` method for easy saving of fit results
   - ``TransmissionModel.load()`` automatically detects model or result files
   - Initialize models directly from saved files: ``TransmissionModel("fit.json")``
   - Loaded results have all methods (plot, save, etc.)
   - JSON-based format avoids ctypes pickle issues
   - See :doc:`save_load` for complete documentation

**Model Parameter Control**
   - New ``vary_basic`` parameter to control thickness and norm fitting
   - Automatic initialization of orientation parameters from material properties
   - Temperature parameter now fixed by default

**MTEX Integration Improvements**
   - ``from_mtex()`` now accepts material as string or dictionary
   - Automatic ``short_name`` generation from material filenames
   - Better handling of powder phases in oriented materials

**Enhanced Documentation**
   - Comprehensive save/load guide
   - Crystal orientation guide
   - Quick reference cards for common operations
   - Executable examples for all major features

Getting Help
------------

- :doc:`/quickstart` - Get started quickly
- :doc:`/api/index` - Complete API reference
- :doc:`/examples/iron_powder` - Working examples
- `GitHub Issues <https://github.com/TsvikiHirsh/nbragg/issues>`_ - Report bugs or request features
