.. models_physics:

Physics of the TransmissionModel
================================

The TransmissionModel in the nres package is based on the fundamental principles of neutron transmission through materials. This model incorporates cross-section calculations, instrument response functions, and background parameterization to accurately describe experimental data.

Transmission Equation
---------------------

The core equation for the transmission of neutrons through a material is given by:

.. math::

   T(E) = \text{norm} \cdot e^{-\sigma(E) \cdot \text{thickness} \cdot n} \cdot (1 - \text{bg}(E)) + \text{bg}(E)

Where:
   - :math:`T(E)` is the transmission as a function of energy
   - :math:`\text{norm}` is a normalization factor
   - :math:`\sigma(E)` is the energy-dependent cross-section
   - :math:`\text{thickness}` is the sample thickness
   - :math:`n` is the total atomic weight (in atoms/barn-cm)
   - :math:`\text{bg}(E)` is the energy-dependent background function

Cross-Section Integration
-------------------------

The cross-section :math:`\sigma(E)` is typically provided as a set of discrete points. To apply this to experimental data, we need to integrate the cross-section over the requested energy bins. This integration is performed using the `integrate_cross_section` function, which employs the following steps:

1. Extend the energy grid to accommodate the convolution kernel.
2. Perform trapezoidal integration over each energy bin.
3. Calculate the average cross-section for each bin.
4. If a response function is provided, convolve the result with the kernel.
5. Trim the result to match the original energy grid size.

This integration ensures that the cross-section is properly averaged over each energy bin, accounting for the energy resolution of the instrument.

Instrument Response Function
----------------------------

The instrument response function accounts for the finite resolution of the measurement apparatus. In the nres package, this is implemented in the `Response` class, which offers two main options:

1. Exponential-Gaussian Response ("expo_gauss"):
   This response function combines an exponential decay with a Gaussian distribution, suitable for many neutron time-of-flight instruments. It is defined as:

   .. math::

      R(t) = \text{exponnorm.pdf}(t, K, \text{loc}=x_0, \text{scale}=\tau)

   Where :math:`K` is the shape parameter, :math:`x_0` is the location parameter, and :math:`\tau` is the scale parameter.

2. Empty Response ("none"):
   This option applies no response function, useful for ideal or theoretical calculations.

The response function is convolved with the integrated cross-section to simulate the instrument's effect on the measured transmission.

Background Parameterization
---------------------------

The background in neutron transmission experiments can arise from various sources, including sample-independent neutron scattering and electronic noise. The `Background` class in nres provides several options for modeling this background:

1. Constant Background:
   .. math::

      \text{bg}(E) = b_0

2. Third-degree Polynomial Background:
   .. math::

      \text{bg}(E) = b_0 + b_1\sqrt{E} + \frac{b_2}{\sqrt{E}}

3. Fifth-degree Polynomial Background:
   .. math::

      \text{bg}(E) = b_0 + b_1\sqrt{E} + \frac{b_2}{\sqrt{E}} + b_3E + b_4E^2

4. No Background ("none"):
   This option assumes zero background, useful for ideal or corrected data.

The choice of background function depends on the specific experimental setup and the energy range of interest.

Model Integration
-----------------

The TransmissionModel combines these components - cross-section, response function, and background - into a cohesive model:

1. The cross-section is integrated over the energy bins.
2. The integrated cross-section is convolved with the instrument response function.
3. The background is calculated for each energy point.
4. The transmission is computed using the equation provided earlier.

This integrated approach allows for accurate modeling of neutron transmission experiments, accounting for the complex interplay between material properties, instrument characteristics, and experimental conditions.

By adjusting the parameters of each component (cross-section weights, response function parameters, background coefficients, etc.), the model can be fitted to experimental data, providing insights into the sample composition and structure.