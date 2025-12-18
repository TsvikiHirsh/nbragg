====================================
CrossSection Parameter Assignment
====================================

This guide demonstrates the ergonomic parameter assignment feature in ``CrossSection`` that allows users to easily add extinction, SANS, and orientation parameters directly in the constructor.

Overview
--------

Instead of using nested dictionaries to specify material parameters, you can now pass parameters as keyword arguments directly in the ``CrossSection`` constructor. Parameters are automatically assigned to materials based on their position in the arguments.

**Key Benefits:**

* Cleaner, more readable code
* Less nesting and boilerplate
* Parameters assigned in order of appearance
* Fully backward compatible

Single Material with Parameters
--------------------------------

The simplest use case is adding parameters to a single material:

Extinction Parameters
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from nbragg import CrossSection

    # Add extinction parameters to a material
    xs = CrossSection(
        steel='Al_sg225.ncmat',
        ext_l=100,      # Crystallite size (microns)
        ext_Gg=1500,    # Mosaic spread (seconds of arc)
        ext_L=2000      # Strain parameter
    )

This is equivalent to the traditional approach:

.. code-block:: python

    # Traditional approach (still supported)
    xs = CrossSection({
        'steel': {
            'mat': 'Al_sg225.ncmat',
            'ext_l': 100,
            'ext_Gg': 1500,
            'ext_L': 2000
        }
    })

SANS Parameters
~~~~~~~~~~~~~~~

.. code-block:: python

    # Add SANS hard-sphere radius
    xs = CrossSection(
        aluminum='Al_sg225.ncmat',
        sans=50.0  # Hard-sphere radius in Angstroms
    )

Orientation Parameters
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Add orientation parameters
    xs = CrossSection(
        sample='Al_sg225.ncmat',
        theta=45,   # Polar angle (degrees)
        phi=90,     # Azimuthal angle (degrees)
        mos=5       # Mosaicity spread (degrees)
    )

Multiple Materials with Order-Based Assignment
-----------------------------------------------

When you have multiple materials, parameters are assigned based on their order in the keyword arguments. Parameters following a material are assigned to that material:

.. code-block:: python

    # ext_Gg goes to mat1, sans goes to mat2
    xs = CrossSection(
        mat1='Al_sg225.ncmat', ext_Gg=100,
        mat2='Al_sg225.ncmat', sans=20
    )

    # Verify the assignments
    print(xs.materials['mat1']['ext_Gg'])  # Output: 100
    print(xs.materials['mat1']['sans'])    # Output: None
    print(xs.materials['mat2']['ext_Gg'])  # Output: None
    print(xs.materials['mat2']['sans'])    # Output: 20

Multiple Parameters per Material
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can assign multiple parameters to each material by placing them after the material definition:

.. code-block:: python

    xs = CrossSection(
        iron='Al_sg225.ncmat', ext_l=50, ext_Gg=200,      # Iron gets extinction params
        aluminum='Al_sg225.ncmat', theta=30, phi=60        # Aluminum gets orientation params
    )

Combining Different Parameter Types
------------------------------------

You can mix and match different types of parameters for comprehensive material configuration:

.. code-block:: python

    xs = CrossSection(
        sample='Al_sg225.ncmat',
        # Extinction parameters
        ext_l=100,
        ext_Gg=1500,
        # Orientation parameters
        theta=45,
        phi=90,
        mos=5,
        # Lattice parameters
        a=4.05,
        b=4.05,
        c=4.05,
        # Other parameters
        temp=400
    )

Supported Parameters
--------------------

The following parameters can be assigned via keyword arguments:

Extinction Parameters
~~~~~~~~~~~~~~~~~~~~~

* ``ext_l`` - Crystallite size (microns)
* ``ext_Gg`` - Mosaic spread (seconds of arc)
* ``ext_L`` - Strain parameter
* ``ext_method`` - Extinction method
* ``ext_dist`` - Distribution type

SANS Parameters
~~~~~~~~~~~~~~~

* ``sans`` - Hard-sphere radius (Angstroms)

Orientation Parameters
~~~~~~~~~~~~~~~~~~~~~~

* ``theta`` - Polar angle (degrees)
* ``phi`` - Azimuthal angle (degrees)
* ``mos`` - Mosaicity spread (degrees)
* ``dir1`` - First direction vector
* ``dir2`` - Second direction vector
* ``dirtol`` - Direction tolerance

Lattice Parameters
~~~~~~~~~~~~~~~~~~

* ``a`` - Lattice constant a (Angstroms)
* ``b`` - Lattice constant b (Angstroms)
* ``c`` - Lattice constant c (Angstroms)

Other Parameters
~~~~~~~~~~~~~~~~

* ``temp`` - Temperature (Kelvin)
* ``weight`` - Material weight (for multi-phase materials)

Complete Example
----------------

Here's a comprehensive example using multiple materials and various parameters:

.. code-block:: python

    from nbragg import CrossSection, TransmissionModel, Data

    # Create CrossSection with multiple materials and parameters
    xs = CrossSection(
        # First material: Iron with extinction
        iron='Fe_sg229_Iron-alpha.ncmat',
        ext_l=100,
        ext_Gg=1500,

        # Second material: Aluminum with orientation
        aluminum='Al_sg225.ncmat',
        theta=30,
        phi=60,
        mos=3
    )

    # Create model with the configured cross-section
    model = TransmissionModel(
        xs,
        vary_basic=True,
        vary_background=True,
        vary_extinction=True,
        vary_orientation=True
    )

    # Load data and fit
    data = Data.from_transmission("measurement.csv")
    result = model.fit(data)

    # Plot results
    result.plot()

Backward Compatibility
----------------------

This feature is **fully backward compatible**. All existing code using the traditional dictionary-based approach continues to work without any changes:

.. code-block:: python

    # Old style - still works perfectly
    xs = CrossSection({
        'material': {
            'mat': 'Al_sg225.ncmat',
            'ext_l': 100,
            'theta': 45
        }
    })

    # New style - more ergonomic
    xs = CrossSection(
        material='Al_sg225.ncmat',
        ext_l=100,
        theta=45
    )

Both approaches produce identical results and can be mixed in the same codebase.

Tips and Best Practices
------------------------

1. **Use descriptive material names** for clarity:

   .. code-block:: python

       xs = CrossSection(
           substrate='Al_sg225.ncmat', ext_l=100,
           coating='Fe_sg229_Iron-alpha.ncmat', sans=50
       )

2. **Group related parameters** for readability:

   .. code-block:: python

       xs = CrossSection(
           sample='Al_sg225.ncmat',
           # Extinction
           ext_l=100, ext_Gg=1500,
           # Orientation
           theta=45, phi=90,
           # Temperature
           temp=300
       )

3. **Parameter order matters** when using multiple materials - parameters are assigned to the material that precedes them

4. **Use the traditional dictionary approach** when you need more complex configurations or when material parameters come from external sources

Running the Example
-------------------

You can run the complete example code:

.. code-block:: bash

    python docs/source/examples/parameter_assignment.py

This will print detailed output showing how parameters are assigned to different materials.

See Also
--------

* :doc:`../user_guide/cross_section` - Complete CrossSection documentation
* :doc:`../user_guide/extinction` - Extinction modeling details
* :doc:`orientation_examples` - Orientation parameter examples
* :doc:`../api/cross_section` - CrossSection API reference
