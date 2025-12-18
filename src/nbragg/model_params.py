"""
Parameter creation methods for TransmissionModel.

This module contains mixin class with methods for creating different types
of parameters used in neutron transmission fitting.
"""

import lmfit
import numpy as np
import warnings


class ParametersMixin:
    """Mixin class providing parameter creation methods for TransmissionModel."""

    def _make_basic_params(self, vary=True):
        """
        Create basic transmission parameters (thickness and norm).

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default True.

        Returns
        -------
        lmfit.Parameters
            The basic parameters (thickness, norm).
        """
        params = lmfit.Parameters()
        params.add("thickness", value=1., min=0., vary=vary)
        params.add("norm", value=1., min=0., vary=vary)
        return params

    def _make_temperature_params(self):
        """
        Create temperature parameter.

        Returns
        -------
        lmfit.Parameters
            The temperature parameter (always fixed by default at 293.15 K).
        """
        params = lmfit.Parameters()
        params.add("temp", value=293.15, min=0., vary=False)  # Default temperature in Kelvin, always fixed by default
        return params

    def _make_weight_params(self, vary=False):
        """
        Create isotope/phase weight parameters.

        For multi-phase materials, creates parameters that ensure weights sum to 1.
        Uses log-ratio parameterization for numerical stability.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The weight parameters.
        """
        params = lmfit.Parameters()
        weights = np.array([self._materials[phase]["weight"] for phase in self._materials])
        param_names = [phase.replace("-", "") for phase in self._materials]

        N = len(weights)
        if N == 1:
            # Special case: if N=1, the weight is always 1
            params.add(f'{param_names[0]}', value=1., vary=False)
        else:

            last_weight = weights[-1]
            # Add (N-1) free parameters corresponding to the first (N-1) items
            for i, name in enumerate(param_names[:-1]):
                initial_value = weights[i]  # Use weight values
                params.add(f'p{i+1}',value=np.log(weights[i]/last_weight),min=-14,max=14,vary=vary) # limit to 1ppm

            # Define the normalization expression
            normalization_expr = ' + '.join([f'exp(p{i+1})' for i in range(N-1)])

            # Add weights based on the free parameters
            for i, name in enumerate(param_names[:-1]):
                params.add(f'{name}', expr=f'exp(p{i+1}) / (1 + {normalization_expr})')

            # The last weight is 1 minus the sum of the previous weights
            params.add(f'{param_names[-1]}', expr=f'1 / (1 + {normalization_expr})')

        return params

    def _make_lattice_params(self, vary=False):
        """
        Create lattice-parameter ('a','b','c') params for the model.

        Automatically detects crystal symmetry and creates appropriate parameter
        constraints (e.g., a=b=c for cubic, a=b≠c for hexagonal, etc.).

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The lattice-related parameters.
        """
        params = lmfit.Parameters()
        for i, material in enumerate(self._materials):
            # update materials with new lattice parameter
            try:
                info = self.cross_section.phases_data[material].info.structure_info
                a, b, c = info["a"], info["b"], info["c"]

                param_a_name = f"a{i+1}" if len(self._materials)>1 else "a"
                param_b_name = f"b{i+1}" if len(self._materials)>1 else "b"
                param_c_name = f"c{i+1}" if len(self._materials)>1 else "c"

                if np.isclose(a,b,atol=1e-4) and np.isclose(b,c,atol=1e-4):
                    if param_a_name in self.params:
                        self.params[param_a_name].vary = vary
                    else:
                        params.add(param_a_name, value=a, min=0.5, max=10, vary=vary)
                        params.add(param_b_name, value=a, min=0.5, max=10, vary=vary, expr=param_a_name)
                        params.add(param_c_name, value=a, min=0.5, max=10, vary=vary, expr=param_a_name)
                elif np.isclose(a,b,atol=1e-4) and not np.isclose(c,b,atol=1e-4):
                    if param_a_name in self.params:
                        self.params[param_a_name].vary = vary
                        self.params[param_c_name].vary = vary
                    else:
                        params.add(param_a_name, value=a, min=0.5, max=10, vary=vary)
                        params.add(param_b_name, value=a, min=0.5, max=10, vary=vary, expr=param_a_name)
                        params.add(param_c_name, value=c, min=0.5, max=10, vary=vary)
                elif not np.isclose(a,b,atol=1e-4) and not np.isclose(c,b,atol=1e-4):
                    if param_a_name in self.params:
                        self.params[param_a_name].vary = vary
                        self.params[param_b_name].vary = vary
                        self.params[param_c_name].vary = vary
                    else:
                        params.add(param_a_name, value=a, min=0.5, max=10, vary=vary)
                        params.add(param_b_name, value=b, min=0.5, max=10, vary=vary)
                        params.add(param_c_name, value=c, min=0.5, max=10, vary=vary)
            except:
                pass

        return params

    def _make_extinction_params(self, vary=False):
        """
        Create extinction-parameter ('ext_l', 'ext_Gg', 'ext_L') params for the model.

        Requires the CrysExtn NCrystal plugin to be installed.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The extinction-related parameters.
        """
        params = lmfit.Parameters()
        for i, material in enumerate(self._materials):
            # update materials with new lattice parameter
            try:
                info = self.cross_section.extinction[material]


                l, Gg, L = info["l"], info["Gg"], info["L"]

                param_l_name = f"ext_l{i+1}" if len(self._materials)>1 else "ext_l"
                param_Gg_name = f"ext_Gg{i+1}" if len(self._materials)>1 else "ext_Gg"
                param_L_name = f"ext_L{i+1}" if len(self._materials)>1 else "ext_L"


                if param_l_name in self.params:
                    self.params[param_l_name].vary = vary
                    self.params[param_Gg_name].vary = vary
                    self.params[param_L_name].vary = vary
                else:
                    params.add(param_l_name, value=l, min=0., max=10000,vary=vary)
                    params.add(param_Gg_name, value=Gg, min=0., max=10000,vary=vary)
                    params.add(param_L_name, value=L, min=0., max=1000000,vary=vary)
            except KeyError:
                warnings.warn(f"@CRYSEXTN section is not defined for the {material} phase")

        return params

    def _make_sans_params(self, vary=False):
        """
        Create SANS hard-sphere radius parameters for the model.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The SANS-related parameters.
        """
        params = lmfit.Parameters()
        for i, material in enumerate(self._materials):
            # Check if material has sans defined
            sans_value = self._materials[material].get('sans')
            if sans_value is not None:
                param_sans_name = f"sans{i+1}" if len(self._materials) > 1 else "sans"

                if param_sans_name in self.params:
                    self.params[param_sans_name].vary = vary
                else:
                    params.add(param_sans_name, value=sans_value, min=0., max=1000, vary=vary)

        return params

    def _make_orientation_params(self, vary=False):
        """
        Create orientation parameters (θ, ϕ, η) for each phase.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The orientation-related parameters.
        """
        params = lmfit.Parameters()
        materials = self.cross_section.materials
        for phase in self.cross_section.phases:
            # Get orientation values from material dictionary, default to 0
            material = materials.get(phase, {})
            theta_val = material.get('theta', 0.) if material.get('theta') is not None else 0.
            phi_val = material.get('phi', 0.) if material.get('phi') is not None else 0.
            mos_val = material.get('mos', 0.) if material.get('mos') is not None else 0.

            params.add(f"θ_{phase}", value=theta_val, vary=vary)
            params.add(f"ϕ_{phase}", value=phi_val, vary=vary)
            params.add(f"η_{phase}", value=mos_val, min=0., vary=vary)
        return params

    def _make_tof_params(self, vary=False, **kwargs):
        """
        Create time-of-flight correction parameters (L0, t0).

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The TOF-related parameters.
        """
        params = lmfit.Parameters()
        params.add("L0", value=1., min=0., max = 2., vary=vary)
        params.add("t0", value=0., vary=vary)
        return params

    def _tof_correction(self, E, **kwargs):
        """
        Apply time-of-flight correction to energy values.

        Parameters
        ----------
        E : array-like
            Energy values.
        **kwargs
            Keyword arguments including L0 and t0.

        Returns
        -------
        array-like
            Corrected energy values.
        """
        L0 = kwargs.get("L0", self.tof_length)
        t0 = kwargs.get("t0", 0.)
        # Assuming energy correction based on TOF
        return E * (L0 / self.tof_length) + t0
