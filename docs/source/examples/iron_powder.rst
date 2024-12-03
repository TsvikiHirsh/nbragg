============
Iron Powder Example
============

This example demonstrates how to use nbragg to analyze neutron transmission through an iron powder sample.

Dataset Overview
----------------

The iron powder example showcases the analysis of Bragg edges in a polycrystalline iron sample. We'll walk through the entire process from data loading to fitting and interpretation.

Data Preparation
----------------

First, load the transmission data:

.. code-block:: python

    import nbragg
    import matplotlib.pyplot as plt

    # Load transmission data
    data = nbragg.Data.from_transmission("iron_powder.csv")

Cross-Section Configuration
---------------------------

We'll use the NCrystal cross-section for alpha-iron:

.. code-block:: python

    # Define iron cross-section
    xs = nbragg.CrossSection.from_material("Fe_sg229_Iron-alpha.ncmat")

Model Creation and Fitting
--------------------------

Create a transmission model with background and response variations:

.. code-block:: python

    # Create transmission model
    model = nbragg.TransmissionModel(
        xs, 
        vary_background=True, 
        vary_response=True
    )

    # Perform fitting
    result = model.fit(data)

Visualization
-------------

Plot the fitting results:

.. code-block:: python

    # Plot results
    result.plot()
    plt.title("Iron Powder Bragg Edge Fitting")
    plt.show()

Key Observations
----------------

- The iron powder sample shows multiple Bragg edges
- Variations in background and instrumental response are accounted for
- The model provides insights into material structure

Additional Notes
----------------

- Ensure you have the correct NCrystal material file
- Calibrate your instrumental response carefully
- The quality of fitting depends on data resolution