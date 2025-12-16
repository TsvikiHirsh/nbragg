==================
Grouped Data Fitting
==================

Overview
--------

nbragg supports analysis of grouped/gridded data, enabling you to:

- Fit spatially-resolved measurements (e.g., imaging data on a 2D grid)
- Analyze sequential measurements (e.g., scans along a line)
- Process multiple samples or regions of interest simultaneously
- Visualize parameter variations across your dataset

All grouped fitting operations support parallel processing and provide intuitive visualization tools.

Supported Data Structures
--------------------------

Three types of grouped data are supported:

**2D Grids**
   Regular spatial grids (e.g., detector pixels, scan positions)

   - Indexed by ``(row, col)`` tuples
   - Visualized as 2D heatmaps/colormeshes
   - Example: ``result.plot_parameter_map("thickness")`` shows a spatial map

**1D Arrays**
   Linear sequences of measurements

   - Indexed by integers (0, 1, 2, ...)
   - Visualized as line plots
   - Example: parameter evolution along a scan direction

**Named Groups**
   Arbitrary collections with custom identifiers

   - Indexed by descriptive strings ("sample_A", "roi_center", etc.)
   - Visualized as bar charts
   - Example: comparing different sample conditions

Loading Grouped Data
---------------------

Basic Usage
~~~~~~~~~~~

Load data using glob patterns::

    import nbragg

    data = nbragg.Data.from_grouped(
        signal="path/to/signal_*.csv",
        openbeam="path/to/openbeam_*.csv",
        L=10,  # sample-detector distance in meters
        tstep=10e-6  # time step in seconds
    )

The indices are automatically extracted from filenames. For 2D data, use naming like:
- ``signal_x0_y0.csv``, ``signal_x0_y1.csv``, ...
- ``signal_row0_col0.csv``, ``signal_row0_col1.csv``, ...

For 1D data:
- ``signal_0.csv``, ``signal_1.csv``, ...
- ``signal_pixel_0.csv``, ``signal_pixel_1.csv``, ...

Alternative Loading Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**From folders** (loads all CSV files)::

    data = nbragg.Data.from_grouped(
        signal="path/to/signal_folder/",
        openbeam="path/to/openbeam_folder/",
        L=10, tstep=10e-6
    )

**From file lists**::

    signal_files = ["sig1.csv", "sig2.csv", "sig3.csv"]
    openbeam_files = ["ob1.csv", "ob2.csv", "ob3.csv"]

    data = nbragg.Data.from_grouped(
        signal=signal_files,
        openbeam=openbeam_files,
        L=10, tstep=10e-6
    )

Data Attributes
~~~~~~~~~~~~~~~

Grouped data objects have these attributes::

    data.is_grouped       # True for grouped data
    data.indices          # List of string indices: ["(0,0)", "(0,1)", ...]
    data.group_shape      # Shape tuple: (3, 3) for 3x3 grid, (5,) for 5-element array
    data.groups           # Dict mapping indices to dataframes

Fitting Grouped Data
--------------------

Basic Fitting
~~~~~~~~~~~~~

Fitting works identically to single datasets, but processes all groups::

    from nbragg import TransmissionModel, CrossSection, materials

    # Define model
    xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
    model = TransmissionModel(xs, vary_basic=True)

    # Fit all groups (automatically parallelized)
    result = model.fit(
        data,
        n_jobs=4,  # Number of parallel workers
        progress_bar=True,  # Show progress
        wlmin=1.5,
        wlmax=5.0
    )

Parallel Processing
~~~~~~~~~~~~~~~~~~~

Grouped fitting uses true multiprocessing to achieve significant speedup on multi-core systems.

The ``n_jobs`` parameter controls parallelization:

- ``n_jobs=1``: Sequential processing
- ``n_jobs=4``: Use 4 parallel workers (good default)
- ``n_jobs=8``: Use 8 parallel workers
- ``n_jobs=-1``: Use all available CPUs

**Typical speedup** (32 groups, 500 wavelength points):

- 2 workers: ~1.7x faster
- 4 workers: ~2.5x faster
- 8 workers: ~3.3x faster

The ``backend`` parameter selects the parallelization strategy:

- ``backend="loky"`` (default): True multiprocessing with separate processes.
  Each worker reconstructs the model independently, enabling full CPU parallelism.
- ``backend="threading"``: Thread-based parallelism. Limited by Python's GIL
  but lower overhead for very fast fits.
- ``backend="sequential"``: No parallelism. Useful for debugging.

Example with explicit backend::

    # True multiprocessing (default, recommended)
    result = model.fit(data, n_jobs=4, backend="loky")

    # Sequential for debugging
    result = model.fit(data, backend="sequential")

.. note::
   The first batch of fits has initialization overhead (~1-3s total) as each worker
   loads NCrystal. Subsequent fits are fast. For small datasets (<10 groups),
   sequential processing may be faster due to this overhead.

Result Structure
~~~~~~~~~~~~~~~~

The ``fit()`` method returns a ``GroupedFitResult`` object::

    result.results        # Dict mapping indices to individual ModelResult objects
    result.indices        # List of group indices
    result.group_shape    # Original data shape

Accessing Results
-----------------

Individual Group Results
~~~~~~~~~~~~~~~~~~~~~~~~

Access results using flexible indexing::

    # Using tuples (for 2D grids)
    r = result[(1, 2)]

    # Using strings (spaces optional)
    r = result["(1,2)"]    # No spaces
    r = result["(1, 2)"]   # With spaces

    # Using integers (for 1D arrays)
    r = result[5]
    r = result["5"]

    # Using names (for named groups)
    r = result["sample_A"]

Each individual result is a standard lmfit ``ModelResult`` with all normal methods::

    r.params['thickness'].value  # Parameter value
    r.params['thickness'].stderr # Parameter error
    r.success                    # Fit success flag
    r.redchi                     # Reduced chi-square
    r.plot()                     # Plot fit

Summary Statistics
~~~~~~~~~~~~~~~~~~

Get a comprehensive summary of all fits::

    summary_df = result.summary()

This returns a pandas DataFrame with columns:

- ``index``: Group identifier
- ``success``: Fit success flag
- ``redchi``, ``chisqr``: Fit statistics
- ``nfev``, ``nvarys``: Evaluation count and parameter count
- ``<param_name>``: Value for each fitted parameter
- ``<param_name>_err``: Error for each fitted parameter

In Jupyter notebooks, display a formatted HTML report::

    from IPython.display import HTML, display
    display(HTML(result.fit_report()))

Plotting Individual Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plot any individual group's fit::

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    result.plot(index=(1, 1), ax=ax)
    plt.show()

For multi-stage fits, view the stages progression::

    stages_table = result.stages_summary(index=(1, 1))
    print(stages_table)

Plot the total cross-section::

    result.plot_total_xs(index=(1, 1), plot_dspace=True)

Visualizing Parameter Maps
---------------------------

The ``plot_parameter_map()`` method is the primary visualization tool. It automatically
detects the appropriate plot type based on your data structure.

Basic Parameter Maps
~~~~~~~~~~~~~~~~~~~~

Plot any fitted parameter::

    result.plot_parameter_map("thickness")

For 2D grids, this creates a heatmap. For 1D arrays, a line plot. For named groups, a bar chart.

Customizing Appearance
~~~~~~~~~~~~~~~~~~~~~~

Control the visualization with keyword arguments::

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))
    result.plot_parameter_map(
        "thickness",
        ax=ax,
        cmap="viridis",        # Colormap (for 2D)
        plot_errors=False,     # Show/hide error map
        vmin=0.9,              # Color scale minimum
        vmax=1.1,              # Color scale maximum
    )
    plt.title("Sample thickness variation")
    plt.tight_layout()
    plt.show()

Plotting Errors
~~~~~~~~~~~~~~~

Visualize parameter uncertainties::

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Values
    result.plot_parameter_map("thickness", ax=ax1)
    ax1.set_title("Thickness values")

    # Errors
    result.plot_parameter_map("thickness", plot_errors=True, ax=ax2)
    ax2.set_title("Thickness errors")

    plt.show()

Filtering with Queries
~~~~~~~~~~~~~~~~~~~~~~

Use pandas query syntax to filter which groups to display::

    # Show only successful fits
    result.plot_parameter_map(
        "thickness",
        query="success == True"
    )

    # Show only good fits (reduced chi-square < 2)
    result.plot_parameter_map(
        "thickness",
        query="success == True and redchi < 2.0"
    )

    # Complex queries
    result.plot_parameter_map(
        "thickness",
        query="redchi < 1.5 and thickness > 0.8 and thickness < 1.2"
    )

Explicit Plot Type Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Override automatic detection if needed::

    # Force 2D heatmap (if you have 2D data but want line plot, this won't work well)
    result.plot_parameter_map("thickness", kind="pcolormesh")

    # Force 1D line plot
    result.plot_parameter_map("thickness", kind="line")

    # Force bar chart
    result.plot_parameter_map("thickness", kind="bar")

Saving and Loading
------------------

Saving Results
~~~~~~~~~~~~~~

Save all grouped fit results to a single file::

    result.save("grouped_results.json")

With optional model saving::

    result.save(
        "grouped_results.json",
        model_filename="model.json"
    )

Compact format (faster but less readable)::

    result.save("grouped_results.json", compact=True)

Loading Results
~~~~~~~~~~~~~~~

Load saved grouped results::

    from nbragg.models import GroupedFitResult

    result = GroupedFitResult.load("grouped_results.json")

Load with a model::

    result = GroupedFitResult.load(
        "grouped_results.json",
        model_filename="model.json"
    )

Or pass a model instance::

    result = GroupedFitResult.load(
        "grouped_results.json",
        model=my_model
    )

Combining Datasets
------------------

Add Multiple Measurements
~~~~~~~~~~~~~~~~~~~~~~~~~

Combine data from multiple measurement runs::

    data1 = nbragg.Data.from_grouped(
        signal="run1/signal_*.csv",
        openbeam="run1/ob_*.csv",
        L=10, tstep=10e-6
    )

    data2 = nbragg.Data.from_grouped(
        signal="run2/signal_*.csv",
        openbeam="run2/ob_*.csv",
        L=10, tstep=10e-6
    )

    # Combine (adds counts, properly propagates errors)
    combined = data1 + data2

.. warning::
   Both datasets must have identical indices and group structure. Mismatched indices raise a ``ValueError``.

Advanced Features
-----------------

Multi-Stage Fitting
~~~~~~~~~~~~~~~~~~~

Use Rietveld-type staged refinement::

    # Define stages
    model.stages = {
        'basic': ['norm', 'thickness'],
        'background': 'background',
        'all': 'all'
    }

    # Fit with stages
    result = model.fit(data, stages='all', n_jobs=4)

    # View stage progression for any group
    result.stages_summary(index=(1, 1))

See :doc:`advanced_fitting` for detailed information on multi-stage fitting.

Custom Index Ordering
~~~~~~~~~~~~~~~~~~~~~

Control the order of indices during loading::

    data = nbragg.Data.from_grouped(
        signal=signal_files,
        openbeam=ob_files,
        indices=["center", "edge", "corner"],  # Custom order
        L=10, tstep=10e-6
    )

Working with Subsets
~~~~~~~~~~~~~~~~~~~~

Extract and fit specific groups::

    # Get subset of data (example - would need custom method)
    # subset_indices = ["(0,0)", "(1,1)", "(2,2)"]
    # subset_data = data.subset(subset_indices)

    # Fit only the subset
    # subset_result = model.fit(subset_data)

Best Practices
--------------

File Naming Conventions
~~~~~~~~~~~~~~~~~~~~~~~

For automatic index extraction:

**2D grids** - Use clear row/column notation::

    signal_x0_y0.csv, signal_x0_y1.csv, signal_x0_y2.csv
    signal_x1_y0.csv, signal_x1_y1.csv, signal_x1_y2.csv

or::

    signal_row0_col0.csv, signal_row0_col1.csv, ...

**1D arrays** - Use sequential numbering::

    signal_0.csv, signal_1.csv, signal_2.csv, ...

or::

    signal_pixel_0.csv, signal_pixel_1.csv, ...

**Named groups** - Use descriptive names::

    signal_sample_A.csv, signal_sample_B.csv
    signal_center.csv, signal_edge.csv, signal_corner.csv

Performance Tips
~~~~~~~~~~~~~~~~

1. **Use appropriate n_jobs**: Start with ``n_jobs=4`` for most systems, increase for large datasets
2. **Consider dataset size**: For <10 groups, ``backend="sequential"`` may be faster due to initialization overhead
3. **Limit wavelength range**: Use ``wlmin`` and ``wlmax`` to focus on relevant regions
4. **Pre-filter failed fits**: Use queries in visualizations to hide bad fits
5. **Save intermediate results**: Save after fitting to avoid recomputation
6. **Monitor speedup**: For best efficiency, ensure each fit takes >100ms; faster fits have more overhead

Memory Management
~~~~~~~~~~~~~~~~~

For very large datasets:

- Fit in batches if memory is limited
- Use ``compact=True`` when saving to reduce file size
- Clear unwanted results from memory::

    del result  # Free memory

Common Patterns
---------------

Imaging Analysis
~~~~~~~~~~~~~~~~

Typical workflow for spatially-resolved imaging data::

    import nbragg
    import matplotlib.pyplot as plt

    # Load 2D grid data
    data = nbragg.Data.from_grouped(
        signal="imaging/signal_x*_y*.csv",
        openbeam="imaging/ob_x*_y*.csv",
        L=10, tstep=10e-6
    )

    # Setup model
    xs = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha"])
    model = nbragg.TransmissionModel(xs, vary_basic=True)

    # Fit all positions
    result = model.fit(data, n_jobs=4, progress_bar=True)

    # Create parameter maps
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    result.plot_parameter_map("thickness", ax=ax1)
    result.plot_parameter_map("thickness", plot_errors=True, ax=ax2)
    result.plot_parameter_map("norm", ax=ax3)
    result.plot_parameter_map("norm", plot_errors=True, ax=ax4)

    plt.tight_layout()
    plt.show()

    # Save results
    result.save("imaging_results.json")

Line Scan Analysis
~~~~~~~~~~~~~~~~~~

For measurements along a single direction::

    # Load 1D array data
    data = nbragg.Data.from_grouped(
        signal="scan/position_*.csv",
        openbeam="scan/ob_*.csv",
        L=10, tstep=10e-6
    )

    # Fit
    result = model.fit(data, n_jobs=4)

    # Plot parameter evolution
    fig, ax = plt.subplots(figsize=(10, 6))
    result.plot_parameter_map("thickness", ax=ax)
    ax.set_xlabel("Position along scan")
    ax.set_ylabel("Thickness (cm)")
    plt.show()

Multi-Sample Comparison
~~~~~~~~~~~~~~~~~~~~~~~

Comparing different samples or conditions::

    # Load named groups
    data = nbragg.Data.from_grouped(
        signal="samples/*_signal.csv",
        openbeam="samples/*_ob.csv",
        L=10, tstep=10e-6
    )

    # Fit
    result = model.fit(data, n_jobs=-1)

    # Compare parameters
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    result.plot_parameter_map("thickness", ax=ax1)
    ax1.set_title("Thickness comparison")

    result.plot_parameter_map("norm", ax=ax2)
    ax2.set_title("Normalization comparison")

    plt.show()

    # Get summary table
    summary = result.summary()
    print(summary[['index', 'thickness', 'thickness_err', 'redchi']])

See Also
--------

- :doc:`basic_usage` - Fundamental concepts
- :doc:`advanced_fitting` - Multi-stage fitting strategies
- :doc:`save_load` - Saving and loading results
- :doc:`model_parameters` - Parameter control

For a complete worked example, see the `grouped fits tutorial notebook <https://github.com/TsvikiHirsh/nbragg/blob/master/notebooks/grouped_fits_tutorial.ipynb>`_.
