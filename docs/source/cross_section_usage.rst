Using nres.CrossSection for Material Definition and Combination
===============================================================

This guide demonstrates various ways to define materials and combine them using the ``nres.CrossSection.from_material`` method. We'll also cover how to view cross-section weights, access the cross-section table, and plot the results.

Defining Materials
------------------

Manually
^^^^^^^^

You can create a ``CrossSection`` object by providing isotope names and associated weight fractions (weights) in a dictionary format:

.. code-block:: python

    import nres
    # Create a CrossSection for Carbon, Hydrogen and Iron isotopic mixture
    xs = nres.CrossSection({"C12":0.3,"H1":0.1,"Fe56":0.6})
    
    # View the weights of the mixture
    print(xs.weights)

Output:

.. code-block:: text

    C12     0.3
    H1      0.1
    Fe56    0.6
    dtype: float64

    # Plot the cross-section
    xs.plot(loglog=True)

Working with nres Dictionaries
------------------------------

You can access predefined materials, elements, and isotopes from nres dictionaries:

.. code-block:: python

    # Using nres.materials
    steel_xs = nres.CrossSection.from_material(nres.materials["Steel, Stainless 304"])
    print(steel_xs.weights)

Output:

.. code-block:: text

    C     0.00080
    Mn    0.02000
    P     0.00045
    S     0.00030
    Si    0.01000
    Cr    0.19000
    Ni    0.09500
    Fe    0.68345
    dtype: float64

    # Using nres.elements
    carbon_xs = nres.CrossSection.from_material(nres.elements["Carbon"])
    print(carbon_xs.weights)

Output:

.. code-block:: text

    C-12    0.9893
    C-13    0.0107
    dtype: float64

    # Using nres.isotopes
    u235_xs = nres.CrossSection.from_material(nres.isotopes["U235"])
    print(u235_xs.weights)

Output:

.. code-block:: text

    U-235    1.0
    dtype: float64

Combining Materials
-------------------

You can combine different materials using the ``__add__`` method:

.. code-block:: python

    # Combine iron and nickel
    iron_xs = nres.CrossSection.from_material("Fe")
    iron_nickel_xs = iron_xs + nres.CrossSection.from_material("Ni")
    
    # View the weights of the combined material
    print(iron_nickel_xs.weights)

Output:

.. code-block:: text

    Fe-54    0.029225
    Fe-56    0.458770
    Fe-57    0.010595
    Fe-58    0.001410
    Ni-58    0.405019
    Ni-60    0.079223
    Ni-61    0.003477
    Ni-62    0.011155
    Ni-64    0.001126
    dtype: float64

    # Plot the combined cross-section
    iron_nickel_xs.plot(title="Iron-Nickel Alloy")

Specifying Split Options
------------------------

The ``from_material`` method allows you to specify how to split the cross-sections:

.. code-block:: python

    # Split by isotopes
    water_isotopes = nres.CrossSection.from_material("H2O", splitby="isotopes")
    print(water_isotopes.weights)

Output:

.. code-block:: text

    H-1     0.111894
    H-2     0.000026
    O-16    0.888002
    O-17    0.000038
    O-18    0.000040
    dtype: float64

    # Split by elements
    water_elements = nres.CrossSection.from_material("H2O", splitby="elements")
    print(water_elements.weights)

Output:

.. code-block:: text

    H    0.111920
    O    0.888080
    dtype: float64

    # Split by materials (useful for complex mixtures)
    water_material = nres.CrossSection.from_material("H2O", splitby="materials")
    print(water_material.weights)

Output:

.. code-block:: text

    H2O    1.0
    dtype: float64

Viewing and Analyzing Cross-Sections
------------------------------------

Accessing Weights and Table Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # View weights of components
    print(steel_xs.weights)

    # Access the cross-section table
    print(steel_xs.table.head())

Output:

.. code-block:: text

                   C         Mn          P          S         Si         Cr         Ni         Fe     total
    energy                                                                                                  
    1.000000e-05  0.003148  0.113892  0.000563  0.000405  0.003176  0.602253  0.470814  3.305265  4.499516
    1.000990e-05  0.003148  0.113870  0.000562  0.000405  0.003176  0.602137  0.470726  3.304605  4.498629
    1.001981e-05  0.003148  0.113847  0.000562  0.000405  0.003175  0.602021  0.470638  3.303945  4.497741
    1.002972e-05  0.003147  0.113825  0.000562  0.000405  0.003175  0.601905  0.470550  3.303286  4.496855
    1.003964e-05  0.003147  0.113802  0.000562  0.000405  0.003174  0.601789  0.470462  3.302627  4.495968

Plotting Cross-Sections
^^^^^^^^^^^^^^^^^^^^^^^

The ``plot`` method allows for customization:

.. code-block:: python

    steel_xs.plot(
        title="Steel Cross-Section",
        xlabel="Energy (eV)",
        ylabel="Cross-Section (barn)",
        lw=2,
        logx=True,
        logy=True
    )

Advanced Usage
--------------

Grouping by Isotopes
^^^^^^^^^^^^^^^^^^^^

You can group cross-sections by isotopes:

.. code-block:: python

    grouped_steel = steel_xs.groupby_isotopes()
    print(grouped_steel.weights)

Output:

.. code-block:: text

    C     0.00080
    Mn    0.02000
    P     0.00045
    S     0.00030
    Si    0.01000
    Cr    0.19000
    Ni    0.09500
    Fe    0.68345
    dtype: float64

    grouped_steel.plot(title="Steel Cross-Section (Grouped by Elements)")

This guide provides a comprehensive overview of using ``nres.CrossSection`` for defining and combining materials, as well as analyzing and visualizing cross-section data. Experiment with different materials and combinations to explore their neutron interaction properties!