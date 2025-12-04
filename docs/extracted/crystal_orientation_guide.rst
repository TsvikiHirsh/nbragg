Crystal Orientation in nBragg
==============================

Overview
--------

Understanding how to specify crystal orientations is fundamental to working with nBragg. This guide explains the coordinate system conventions and how to define crystal orientations for Bragg edge analysis.

Coordinate System Convention
-----------------------------

nBragg uses a fixed right-handed orthogonal laboratory coordinate system (:math:`\mathbf{x}`, :math:`\mathbf{y}`, :math:`\mathbf{z}`) where the incident neutron beam propagation vector :math:`\mathbf{k}_i` is always directed along the positive :math:`z`-axis. This simplification makes orientation specification more intuitive compared to the fully flexible coordinate system available in NCrystal.

Defining Crystal Orientation
-----------------------------

Basic Concept
~~~~~~~~~~~~~

The crystal orientation is defined by specifying which reciprocal lattice vectors align with the laboratory coordinate axes. Specifically, nBragg requires you to specify:

* **[hkl]_z**: The Miller indices of the reciprocal lattice vector pointing along the positive :math:`z`-direction (beam direction)
* **[hkl]_y**: The Miller indices of the reciprocal lattice vector pointing along the positive :math:`y`-direction

where :math:`h`, :math:`k`, :math:`l` are Miller indices in reciprocal lattice space.

Mathematical Framework
~~~~~~~~~~~~~~~~~~~~~~

Let :math:`(\mathbf{a}, \mathbf{b}, \mathbf{c})` be the right-handed orthogonal crystal coordinate system. The orientation is described by a rotation :math:`\mathbf{g} \in \text{SO(3)}` that rotates the laboratory coordinate system onto the crystal coordinate system:

.. math::

   g(\alpha, \beta, \gamma) = R(\mathbf{z}, \alpha) R(\mathbf{y}, \beta) R(\mathbf{z}, \gamma)

This is written in Matthies convention, where :math:`\alpha, \gamma \in [0, 2\pi]` and :math:`\beta \in [0, \pi]` are Euler angles.

The unit vector :math:`\hat{\mathbf{n}}_{hkl}` normal to a crystal plane is defined as:

.. math::

   \hat{\mathbf{n}}_{hkl} = h\hat{\mathbf{a}}_1 + k\hat{\mathbf{a}}_2 + l\hat{\mathbf{a}}_3

where :math:`(\hat{\mathbf{a}}_1, \hat{\mathbf{a}}_2, \hat{\mathbf{a}}_3)` are basis vectors of the crystal lattice structure in reciprocal space.

In terms of Euler angles:

.. math::

   \hat{\mathbf{n}}_{hkl} = R(\alpha, \beta, \gamma) \hat{\mathbf{a}}

Practical Examples
------------------

Example 1: Simple Cubic Crystal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a simple cubic crystal with the (001) plane perpendicular to the beam:

.. code-block:: python

   # [hkl]_z = [0, 0, 1] - (001) plane normal to beam
   # [hkl]_y = [0, 1, 0] - (010) plane normal to y-axis
   
   orientation_z = [0, 0, 1]
   orientation_y = [0, 1, 0]

Example 2: Rotated Crystal
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a crystal rotated 45° around the z-axis:

.. code-block:: python

   # [hkl]_z = [0, 0, 1] - still along beam direction
   # [hkl]_y = [1, 1, 0] - now at 45° in xy-plane
   
   orientation_z = [0, 0, 1]
   orientation_y = [1, 1, 0]

Mosaicity and Crystal Imperfections
------------------------------------

Real crystals are imperfect and characterized by their **mosaicity** - the angular spread of crystal orientations around a nominal normal vector :math:`\hat{n}`.

NCrystal models this as a Gaussian distribution with:

* **FWHM** (Full Width at Half Maximum): :math:`\eta`
* **Truncation cutoff**: :math:`\tau = \max(3, 1.1\sqrt{-2\log_e \epsilon}) \eta / (2\sqrt{2\log_e 2})`

where :math:`1 - \epsilon` is the volume fraction within :math:`\tau`, with a default value of :math:`\epsilon = 10^{-3}` in NCrystal.

Crystal Rotation for Alignment
-------------------------------

Understanding Bragg Circles
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Bragg condition is satisfied when a Bragg circle (defined by the neutron wavelength) intersects with the normal vector :math:`\hat{n}` of a crystal plane. For a mosaic crystal, this occurs when the Bragg circle intersects with the distribution of orientations.

Key observations:

* At longer wavelengths (:math:`\lambda_1 > \lambda_2`), Bragg circles translate along the :math:`z`-axis
* **Bragg dips** (single crystals) and **Bragg edges** (polycrystals at :math:`\lambda = 2d_{hkl}`) appear in transmission curves
* Unlike Bragg edges, **Bragg dip positions depend on crystal orientation**

Rotation Parameters
~~~~~~~~~~~~~~~~~~~

nBragg allows you to rotate the crystal around the :math:`x` and :math:`y` axes to align modeled Bragg dips with experimental data:

* :math:`\phi`: Rotation around :math:`x`-axis
* :math:`\theta`: Rotation around :math:`y`-axis

The rotated normal vector is calculated as:

.. math::

   \hat{\mathbf{n}'} = R_Z(\alpha) R_Y(\beta) R_Z(\gamma) R_Y(\theta) R_X(\phi) \hat{\mathbf{a}}

In matrix form:

.. math::

   \begin{bmatrix} h' \\ k' \\ l' \end{bmatrix}
   =
   \begin{bmatrix}
   \cos\theta & 0 & \sin\theta \\
   0 & 1 & 0 \\
   -\sin\theta & 0 & \cos\theta
   \end{bmatrix}
   \begin{bmatrix}
   \cos\phi & -\sin\phi & 0 \\
   \sin\phi & \cos\phi & 0 \\
   0 & 0 & 1
   \end{bmatrix}
   \begin{bmatrix} h \\ k \\ l \end{bmatrix}

Usage in nBragg
~~~~~~~~~~~~~~~

.. code-block:: python

   # Example: Rotating crystal orientation for fitting
   import nbragg
   
   # Initial orientation
   hkl_z = [0, 0, 1]
   hkl_y = [0, 1, 0]
   
   # Rotation angles (in radians or degrees, check your API)
   phi = 0.1    # Rotation around x-axis
   theta = 0.2  # Rotation around y-axis
   
   # Create or update crystal with orientation
   # (Exact API may vary - consult nbragg API documentation)

Visualization and Interpretation
---------------------------------

The coordinate system can be visualized as follows:

* **Incident beam**: Along positive :math:`z`-axis (:math:`\hat{k}_i`)
* **Laboratory frame**: Right-handed :math:`(x, y, z)` system
* **Crystal orientation**: Defined by :math:`[hkl]_z` and :math:`[hkl]_y`
* **Rotated orientation**: Obtained through :math:`\phi` and :math:`\theta` rotations, yielding :math:`[h', k', l']_z` and :math:`[h', k', l']_y`
* **Bragg circles**: Represented as rings at specific wavelengths (:math:`\lambda_1`, :math:`\lambda_2`)
* **xz-plane**: Often highlighted as a reference plane

Key Points for Users
---------------------

Important Conventions
~~~~~~~~~~~~~~~~~~~~~

1. **Fixed beam direction**: Unlike NCrystal's flexibility, nBragg fixes :math:`\mathbf{k}_i` along :math:`+z`
2. **Two-vector specification**: Always specify both :math:`[hkl]_z` and :math:`[hkl]_y`
3. **Right-handed system**: Ensure your vectors maintain right-handedness
4. **Reciprocal space**: Miller indices refer to reciprocal lattice vectors

Common Pitfalls
~~~~~~~~~~~~~~~

* **Orthogonality**: Ensure :math:`[hkl]_z` and :math:`[hkl]_y` are orthogonal in reciprocal space
* **Normalization**: Vectors are automatically normalized, but be aware of this
* **Sign conventions**: Pay attention to positive directions
* **Units**: Check whether rotation angles are in radians or degrees

Best Practices
~~~~~~~~~~~~~~

1. Start with simple orientations (e.g., low-index planes along axes)
2. Verify orientation visually if possible
3. Use small rotation increments when fitting
4. Consider crystal symmetry when choosing orientation
5. Account for mosaicity in your models

Further Reading
---------------

* Matthies, S. (1987). On the diffracting volumes used in pole-figure analysis.
* Schmitt, B. et al. (2023). NCrystal documentation on crystal orientations.
* Kittelmann, T. et al. (2021). NCrystal paper on mosaicity modeling.
* Cai, W. et al. (2020). Bragg-edge neutron transmission imaging.

See Also
--------

* :doc:`nbragg_tutorial` - Basic tutorial for nBragg
* :doc:`api_reference` - Complete API documentation
* :doc:`examples` - Working examples with different crystal systems
