.. models_usage:

TransmissionModel
=================

The ``TransmissionModel`` class in the ``nres`` package is used to define and fit transmission models for nuclear resonance data. It combines cross-section calculations with response functions and background models to provide a comprehensive tool for analyzing transmission spectra.

Basic Usage
-----------

To use the ``TransmissionModel``, you first need to import the necessary modules and create a cross-section object:

.. code-block:: python

   import numpy as np
   from nres.cross_section import CrossSection
   from nres.models import TransmissionModel

   # Create a cross-section object (example)
   xs = CrossSection("Ta181", T=300)

   # Create a TransmissionModel
   model = TransmissionModel(xs, response="expo_gauss", background="polynomial3")

Changing Parameters
-------------------

You can modify the model parameters using the ``model.params`` attribute:

.. code-block:: python

   # Change the thickness parameter
   model.params['thickness'].value = 0.1
   model.params['thickness'].min = 0.01
   model.params['thickness'].max = 1.0

   # Adjust response function parameters
   model.params['sigma'].value = 1e-5
   model.params['tau'].value = 1e-6

   # Modify background parameters
   model.params['b0'].value = 0.01
   model.params['b1'].value = 1e-6

Conducting a Fit
----------------

To fit the model to experimental data, you can use the ``fit`` method:

.. code-block:: python

   import pandas as pd

   # Load your experimental data (example)
   data = pd.read_csv("experimental_data.csv")

   # Perform the fit
   result = model.fit(data, emin=0.5e6, emax=20e6)

   # Print fit report
   print(result.fit_report())

Showing Results and Plotting
----------------------------

After fitting, you can visualize the results using the built-in plotting method:

.. code-block:: python

   # Plot the fit results
   fig, axes = result.plot(plot_bg=True)

   # Customize the plot (optional)
   axes[0].set_ylim(0, 1)
   axes[0].set_xscale('log')

   # Display the plot
   import matplotlib.pyplot as plt
   plt.show()

.. image:: transmission_fit_plot.png
   :alt: Transmission Model Fit Plot
   :width: 600px
   :align: center

The plot shows the experimental data points, the best-fit curve, and the residuals. If ``plot_bg=True``, it also displays the background function.

Advanced Usage
--------------

Varying Weights and Time-of-Flight Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can allow certain parameters to vary during fitting:

.. code-block:: python

   model = TransmissionModel(xs, 
                             vary_weights=True, 
                             vary_tof=True, 
                             vary_background=True, 
                             vary_response=True)

   # Fit with varying parameters
   result = model.fit(data, emin=0.5e6, emax=20e6)

Updating Cross-Section
^^^^^^^^^^^^^^^^^^^^^^

You can update the cross-section of an existing model:

.. code-block:: python

   new_xs = CrossSection("W182", T=300)
   updated_model = model.set_cross_section(new_xs, inplace=False)

Varying All Parameters
^^^^^^^^^^^^^^^^^^^^^^

To quickly vary or fix all parameters:

.. code-block:: python

   # Vary all parameters except 'thickness'
   model.vary_all(vary=True, except_for=['thickness'])

   # Fix all parameters
   model.vary_all(vary=False)

This documentation provides an overview of the ``TransmissionModel`` class and its key functionalities. Users can refer to this guide to understand how to create models, adjust parameters, perform fits, and visualize results using your ``nres`` package.