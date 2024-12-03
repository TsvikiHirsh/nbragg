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
    xs = nbragg.CrossSection.from_material(nbragg.materials["Silicon"])

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

1. **Data Loading**: Use ``nbragg.Data.from_transmission()`` to load experimental data.
2. **Cross-Section**: Define material properties using ``nbragg.CrossSection``.
3. **Model Creation**: Build a transmission model with flexible parameters.
4. **Fitting**: Use the ``fit()`` method to analyze your data.
5. **Visualization**: Easily plot results with the ``plot()`` method.