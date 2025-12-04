Crystal Orientation Quick Reference
====================================

Coordinate System
-----------------

**Laboratory Frame (Fixed)**
  * Right-handed orthogonal system: :math:`(x, y, z)`
  * Incident neutron beam: :math:`\mathbf{k}_i` along :math:`+z` direction
  * **Key difference from NCrystal**: Beam direction is FIXED in nBragg

Specifying Orientation
-----------------------

**Two-Vector Specification**

You must specify TWO reciprocal lattice vectors:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Parameter
     - Description
     - Example
   * - **[hkl]_z**
     - Miller indices along :math:`+z` (beam direction)
     - ``[0, 0, 1]``
   * - **[hkl]_y**
     - Miller indices along :math:`+y` direction
     - ``[0, 1, 0]``

**Requirements**
  * Vectors must be orthogonal in reciprocal space
  * Use Miller indices :math:`(h, k, l)` for reciprocal lattice
  * Right-handed system must be maintained

Common Orientations
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - Description
     - [hkl]_z
     - [hkl]_y
     - Notes
   * - Standard cubic
     - [0, 0, 1]
     - [0, 1, 0]
     - Most common
   * - 45° rotation (z-axis)
     - [0, 0, 1]
     - [1, 1, 0]
     - Rotated in xy-plane
   * - 90° rotation (z-axis)
     - [0, 0, 1]
     - [1, 0, 0]
     - x-axis along y
   * - Arbitrary
     - [1, 1, 2]
     - [1, -1, 0]
     - Check orthogonality

Rotation for Alignment
----------------------

**Rotation Parameters**

When fitting to experimental data:

.. math::

   \text{Rotate around } x\text{-axis: } \phi
   
   \text{Rotate around } y\text{-axis: } \theta

**Rotation Matrix Application**

.. math::

   \begin{bmatrix} h' \\ k' \\ l' \end{bmatrix}
   = R_Y(\theta) \cdot R_X(\phi) \cdot
   \begin{bmatrix} h \\ k \\ l \end{bmatrix}

Where:

.. math::

   R_Y(\theta) = 
   \begin{bmatrix}
   \cos\theta & 0 & \sin\theta \\
   0 & 1 & 0 \\
   -\sin\theta & 0 & \cos\theta
   \end{bmatrix}

.. math::

   R_X(\phi) = 
   \begin{bmatrix}
   1 & 0 & 0 \\
   0 & \cos\phi & -\sin\phi \\
   0 & \sin\phi & \cos\phi
   \end{bmatrix}

Mosaicity Parameters
--------------------

**Gaussian Distribution Model**

.. list-table::
   :widths: 30 70

   * - **η (FWHM)**
     - Full Width at Half Maximum of orientation spread
   * - **τ (cutoff)**
     - :math:`\max(3, 1.1\sqrt{-2\log_e \epsilon}) \cdot \eta / (2\sqrt{2\log_e 2})`
   * - **ε (epsilon)**
     - Volume fraction parameter (default: :math:`10^{-3}`)

**Physical Interpretation**
  * η defines the angular spread of crystal planes
  * Larger η = more mosaic (imperfect) crystal
  * τ determines where distribution is truncated

Bragg Condition
---------------

**Bragg Circles**
  * At wavelength λ, forms a circle around z-axis
  * Longer wavelength → circle moves along z
  * Bragg condition satisfied when circle intersects :math:`\hat{n}`

**For Single Crystals**
  * Sharp Bragg **dips** in transmission
  * Position depends on orientation
  * Use rotation to align dips with data

**For Polycrystals**
  * Bragg **edges** at :math:`\lambda = 2d_{hkl}`
  * Position independent of orientation
  * Edges are orientation-averaged

Troubleshooting Checklist
--------------------------

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Issue
     - Solution
   * - Vectors not orthogonal
     - Check dot product in reciprocal space
   * - Wrong handedness
     - Verify :math:`\mathbf{z} = \mathbf{x} \times \mathbf{y}`
   * - Rotation not working
     - Check angle units (radians vs degrees)
   * - Unexpected Bragg dips
     - Verify crystal structure and orientation
   * - Poor fit to data
     - Adjust mosaicity parameters η

Quick Equations Reference
--------------------------

**Normal Vector**

.. math::

   \hat{\mathbf{n}}_{hkl} = h\hat{\mathbf{a}}_1 + k\hat{\mathbf{a}}_2 + l\hat{\mathbf{a}}_3

**With Euler Angles**

.. math::

   \hat{\mathbf{n}}_{hkl} = R_Z(\alpha) R_Y(\beta) R_Z(\gamma) \hat{\mathbf{a}}

**Rotation (Matthies Convention)**

.. math::

   g(\alpha, \beta, \gamma) = R_Z(\alpha) R_Y(\beta) R_Z(\gamma)

Where: :math:`\alpha, \gamma \in [0, 2\pi]`, :math:`\beta \in [0, \pi]`

Best Practices
--------------

1. **Start Simple**
   
   * Begin with low Miller indices (e.g., [0,0,1], [0,1,0])
   * Test with known crystal orientations

2. **Verify Orthogonality**
   
   * Always check that your chosen vectors are orthogonal
   * Use reciprocal lattice metric for non-cubic systems

3. **Small Increments**
   
   * Use small rotation angles (1-5°) when fitting
   * Avoid large jumps in orientation space

4. **Consider Symmetry**
   
   * Account for crystal symmetry in your analysis
   * Equivalent orientations may exist

5. **Mosaicity Matters**
   
   * Don't ignore mosaicity - it affects results
   * Typical values: 0.1° to 2° for most crystals

Code Snippets
-------------

**Basic Setup**

.. code-block:: python

   # Define orientation
   hkl_z = [0, 0, 1]  # Beam direction
   hkl_y = [0, 1, 0]  # Y direction

**Check Orthogonality**

.. code-block:: python

   import numpy as np
   
   def check_orthogonal(hkl1, hkl2):
       v1 = np.array(hkl1)
       v2 = np.array(hkl2)
       v1 = v1 / np.linalg.norm(v1)
       v2 = v2 / np.linalg.norm(v2)
       return abs(np.dot(v1, v2)) < 1e-10

**Apply Rotation**

.. code-block:: python

   def rotate_orientation(hkl, phi, theta):
       """Rotate [hkl] by angles phi (x) and theta (y)"""
       Rx = np.array([[1, 0, 0],
                      [0, np.cos(phi), -np.sin(phi)],
                      [0, np.sin(phi), np.cos(phi)]])
       
       Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                      [0, 1, 0],
                      [-np.sin(theta), 0, np.cos(theta)]])
       
       return Ry @ Rx @ np.array(hkl)

Related Resources
-----------------

* **Full Documentation**: :doc:`crystal_orientation_guide`
* **Examples**: :doc:`orientation_examples`
* **API Reference**: :doc:`api_reference`
* **Tutorial**: :doc:`nbragg_tutorial`

Citation
--------

When using nBragg's orientation specification in publications:

* Matthies (1987) - Euler angle conventions
* NCrystal documentation - Mosaicity modeling
* Kittelmann et al. (2021) - NCrystal implementation
