==========
Quickstart
==========

This guide will help you get started with nbragg quickly.

Basic Usage
-----------

Here's a simple example of how to use nbragg:

.. code-block:: python

    import nbragg

    # Read transmission data
    data = nbragg.Data.from_transmission("your_data.csv")

    # Define cross-section
    xs = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha.ncmat"])

    # Create transmission model
    model = nbragg.TransmissionModel(xs, 
                                     vary_background=True, 
                                     vary_response=True)

    # Perform fitting
    result = model.fit(data)

    # Plot results
    result.plot()

Key Concepts
------------

1. **Data Loading**: Use ``nbragg.Data.from_transmission()`` to load experimental data
   (or ``nbragg.Data.from_counts()`` to build transmission from sample and open-beam counts).
2. **Cross-Section**: Define material properties using ``nbragg.CrossSection``.
   Combine phases with simple arithmetic: ``xs = 0.5 * alpha + 0.5 * gamma``.
3. **Model Creation**: Build a ``TransmissionModel`` with flexible ``vary_*`` switches
   for background, response, weights, orientation, and more.
4. **Fitting**: Use the ``fit()`` method to analyze your data, optionally with
   multi-stage (Rietveld-style) refinement.
5. **Visualization**: Easily plot results with the ``plot()`` method.

Next Steps
----------

- :doc:`user_guide/basic_usage` — the fundamentals in more depth
- :doc:`user_guide/grouped_fitting` — spatially-resolved and multi-sample data
- :doc:`user_guide/orientation_index` — oriented and textured materials
- :doc:`examples/iron_powder` — a complete worked example