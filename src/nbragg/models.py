import lmfit
import numpy as np
import nbragg.utils as utils
from nbragg.response import Response, Background
from scipy.ndimage import convolve1d
from nbragg.cross_section import CrossSection
from nbragg.data import Data
import NCrystal as NC
import pandas
import matplotlib.pyplot as plt
from copy import deepcopy
from typing import List, Optional, Union, Dict
import warnings
import fnmatch
import re
from numpy import log
import json
import os

# Import refactored components
from nbragg.grouped_fit import GroupedFitResult, _add_save_method_to_result
from nbragg.model_params import ParametersMixin
from nbragg.model_fitting import FittingMixin
from nbragg.model_plotting import PlottingMixin
from nbragg.model_io import IOMixin, save_result, load_result


class TransmissionModel(ParametersMixin, FittingMixin, PlottingMixin, IOMixin, lmfit.Model):
    def __init__(self, cross_section,
                params: "lmfit.Parameters" = None,
                response: str = "jorgensen",
                background: str = "polynomial3",
                tof_length: float = 9,
                vary_basic: bool = None,
                vary_weights: bool = None,
                vary_background: bool = None,
                vary_tof: bool = None,
                vary_response: bool = None,
                vary_orientation: bool = None,
                vary_lattice: bool = None,
                vary_extinction: bool = None,
                vary_sans: bool = None,
                **kwargs):
        """
        Initialize the TransmissionModel, a subclass of lmfit.Model.

        Parameters
        ----------
        cross_section : callable or str
            A CrossSection object, OR a path to a saved model/result JSON file.
            If a string is provided, the model will be loaded from that file.
        response : str, optional
            The type of response function to use, by default "jorgensen".
        background : str, optional
            The type of background function to use, by default "polynomial3".
        tof_length : float, optional
            The flight path length in [m]
        vary_basic : bool, optional
            If True, allows the basic parameters (thickness, norm) to vary during fitting.
            Note: temp parameter is always set to vary=False by default.
        vary_weights : bool, optional
            If True, allows the isotope weights to vary during fitting.
        vary_background : bool, optional
            If True, allows the background parameters (b0, b1, b2) to vary during fitting.
        vary_tof : bool, optional
            If True, allows the TOF (time-of-flight) parameters (L0, t0) to vary during fitting.
        vary_response : bool, optional
            If True, allows the response parameters to vary during fitting.
        vary_orientation : bool, optional
            If True, allows the orientation parameters (θ,ϕ,η) to vary during fitting.
        vary_lattice: bool, optional
            If True, allows the lattice parameters of the material to be varied
        vary_extinction: bool, optional
            If True, allows the extinction parameters of the material to be varied (requires the CrysExtn plugin to be installed)
        vary_sans: bool, optional
            If True, allows the SANS hard-sphere radius parameter to be varied
        kwargs : dict, optional
            Additional keyword arguments for model and background parameters.

        Notes
        -----
        This model calculates the transmission function as a combination of
        cross-section, response function, and background. The fitting stages are automatically
        populated based on the vary_* parameters.

        Examples
        --------
        >>> # Create from CrossSection
        >>> xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
        >>> model = TransmissionModel(xs, vary_background=True)
        >>>
        >>> # Load from saved file
        >>> model = TransmissionModel("my_model.json")
        >>>
        >>> # Load from saved result file
        >>> model = TransmissionModel("my_result.json")
        >>> model.result.plot()  # Access the loaded result
        """
        # Check if cross_section is a file path
        if isinstance(cross_section, str) and os.path.isfile(cross_section):
            # Load from file - need to bypass normal init
            # We'll set a flag and handle this specially
            loaded = self.__class__.load(cross_section)
            # Copy all attributes from loaded model to this instance
            for key, value in loaded.__dict__.items():
                setattr(self, key, value)
            # Don't continue with normal initialization
            return

        # Normal initialization
        super().__init__(self.transmission, **kwargs)

        # make a new instance of the cross section
        self.cross_section = CrossSection(cross_section,
                                        name=cross_section.name,
                                        total_weight=cross_section.total_weight)
        # update atomic density
        self.cross_section.atomic_density = cross_section.atomic_density                                          
        self._materials = self.cross_section.materials
        self.tof_length = tof_length

        if params is not None:
            self.params = params.copy()
        else:
            self.params = lmfit.Parameters()
        if "thickness" not in self.params and "norm" not in self.params:
            if vary_basic is not None:
                self.params += self._make_basic_params(vary=vary_basic)
            else:
                self.params += self._make_basic_params()
        if "temp" not in self.params:
            self.params += self._make_temperature_params()
        if vary_weights is not None:
            self.params += self._make_weight_params(vary=vary_weights)
        if vary_tof is not None:
            self.params += self._make_tof_params(vary=vary_tof, **kwargs)
        if vary_lattice is not None:
            self.params += self._make_lattice_params(vary=vary_lattice)
        if vary_extinction is not None:
            self.params += self._make_extinction_params(vary=vary_extinction)
        if vary_sans is not None:
            self.params += self._make_sans_params(vary=vary_sans)

        self.response = None
        if vary_response is not None:
            self.response = Response(kind=response, vary=vary_response)
            if list(self.response.params.keys())[0] in self.params:
                for param_name in self.params.keys():
                    self.params[param_name].vary = vary_response 
            else:
                self.params += self.response.params

        self.background = None
        if vary_background is not None:
            self.background = Background(kind=background, vary=vary_background)
            if "b0" in self.params:
                for param_name in self.background.params.keys():
                    self.params[param_name].vary = vary_background 
            else:
                self.params += self.background.params

        self.orientation = None
        if vary_orientation is not None:
            self.params += self._make_orientation_params(vary=vary_orientation)

        # set the total atomic weight n [atoms/barn-cm]
        self.atomic_density = self.cross_section.atomic_density

        # Initialize stages based on vary_* parameters
        self._stages = {}
        possible_stages = [
            "basic", "background", "tof", "lattice",
            "mosaicity", "thetas", "phis", "angles", "orientation", "weights", "response", "extinction", "sans"
        ]
        vary_flags = {
            "basic": True if vary_basic is None else vary_basic,  # Default True for backward compatibility
            "background": vary_background,
            "tof": vary_tof,
            "lattice": vary_lattice,
            "mosaicity": vary_orientation,
            "thetas": vary_orientation,
            "phis": vary_orientation,
            "angles": vary_orientation,
            "orientation": vary_orientation,
            "weights": vary_weights,
            "response": vary_response,
            "extinction": vary_extinction,
            "sans": vary_sans,
        }
        for stage in possible_stages:
            if vary_flags.get(stage, False) is True:
                self._stages[stage] = stage

    def transmission(self, wl: np.ndarray, thickness: float = 1, norm: float = 1., **kwargs):
        """
        Transmission function model with background components.

        Parameters
        ----------
        wl : np.ndarray
            The wavelength values at which to calculate the transmission.
        thickness : float, optional
            The thickness of the material (in cm), by default 1.
        norm : float, optional
            Normalization factor, by default 1.
        kwargs : dict, optional
            Additional arguments for background, response, or cross-section.

        Returns
        -------
        np.ndarray
            The calculated transmission values.

        Notes
        -----
        This function combines the cross-section with the response and background
        models to compute the transmission, which is given by:

        .. math:: T(\\lambda) = \\text{norm} \\cdot e^{- \\sigma \\cdot \\text{thickness} \\cdot n} \\cdot (1 - \\text{bg}) + \\text{bg}
        
        where `sigma` is the cross-section, `bg` is the background function, and `n` is the total atomic weight.
        """
        verbose = kwargs.get("verbose",None)
        if verbose:
            print(kwargs)
        E = NC.wl2ekin(wl)
        E = self._tof_correction(E,**kwargs)
        wl = NC.ekin2wl(E)

        if self.background != None:
            k = kwargs.get("k",1.) # sample dependent background factor (k*B)
            bg = self.background.function(wl,**kwargs)
            
        else:
            k = 1.
            bg = 0.

        n = self.atomic_density

        # Transmission function

        xs = self.cross_section(wl,**kwargs)

        if self.response != None:
            response = self.response.function(**kwargs)
            xs = convolve1d(xs,response,0)

        T = norm * np.exp(- xs * thickness * n) * (1 - bg) + k*bg
        return T
    @property
    def stages(self) -> Dict[str, Union[str, List[str]]]:
        """Get the current fitting stages."""
        return self._stages

    @stages.setter
    def stages(self, value: Union[str, Dict[str, Union[str, List[str]]]]):
        """
        Set the fitting stages.

        Parameters
        ----------
        value : str or dict
            If str, must be "all" to use all vary=True parameters.
            If dict, keys are stage names, values are stage definitions ("all", a valid group name, or a list of parameters/groups).
        """
        # Define valid group names from group_map
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
            "lattice": [p for p in self.params if p in ["a", "b", "c"] or p.startswith("a_") or p.startswith("b_") or p.startswith("c_")],
            "extinction": [p for p in self.params if p.startswith("ext_")],
            "sans": [p for p in self.params if p == "sans" or re.compile(r"sans\d+").match(p) or p.startswith("sans_")],
            "orientation": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ") or p.startswith("η")],
            "mosaicity": [p for p in self.params if p.startswith("η")],
            "thetas": [p for p in self.params if p.startswith("θ")],
            "phis": [p for p in self.params if p.startswith("ϕ")],
            "angles": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ")],
            "temperature": [p for p in ["temp"] if p in self.params],
        }
        
        if isinstance(value, str):
            if value != "all":
                raise ValueError("If stages is a string, it must be 'all'")
            self._stages = {"all": "all"}
        elif isinstance(value, dict):
            # Validate stage definitions
            for stage_name, stage_def in value.items():
                if not isinstance(stage_name, str):
                    raise ValueError(f"Stage names must be strings, got {type(stage_name)}")
                if isinstance(stage_def, str):
                    if stage_def != "all" and stage_def not in group_map:
                        raise ValueError(f"Stage definition for '{stage_name}' must be 'all' or a valid group name, got '{stage_def}'")
                elif isinstance(stage_def, list):
                    for param in stage_def:
                        if not isinstance(param, str):
                            raise ValueError(f"Parameters in stage '{stage_name}' must be strings, got {type(param)}")
                else:
                    raise ValueError(f"Stage definition for '{stage_name}' must be 'all', a valid group name, or a list, got {type(stage_def)}")
            self._stages = value
        else:
            raise ValueError(f"Stages must be a string ('all') or dict, got {type(value)}")
    def _repr_html_(self):
        """HTML representation for Jupyter, including parameters and expanded stages tables."""
        from IPython.display import HTML
        import pandas as pd

        # Parameters table
        param_data = []
        for name, param in self.params.items():
            param_data.append({
                'Parameter': name,
                'Value': f"{param.value:.6g}",
                'Vary': param.vary,
                'Min': f"{param.min:.6g}" if param.min is not None else '-inf',
                'Max': f"{param.max:.6g}" if param.max is not None else 'inf',
                'Expr': param.expr if param.expr else ''
            })
        param_df = pd.DataFrame(param_data)
        param_html = param_df.to_html(index=False, classes='table table-striped', border=0)

        # Helper function to resolve a single parameter or group
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
            "lattice": [p for p in self.params if p in ["a", "b", "c"] or p.startswith("a_") or p.startswith("b_") or p.startswith("c_")],
            "extinction": [p for p in self.params if p.startswith("ext_")],
            "sans": [p for p in self.params if p == "sans" or re.compile(r"sans\d+").match(p) or p.startswith("sans_")],
            "orientation": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ") or p.startswith("η")],
            "mosaicity": [p for p in self.params if p.startswith("η")],
            "thetas": [p for p in self.params if p.startswith("θ")],
            "phis": [p for p in self.params if p.startswith("ϕ")],
            "angles": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ")],
            "temperature": [p for p in ["temp"] if p in self.params],
        }

        def resolve_single_param_or_group(item):
            if item == "all":
                return [p for p in self.params if self.params[p].vary]
            elif item in group_map:
                return group_map[item]
            elif item in self.params:
                return [item]
            else:
                matching_params = [p for p in self.params.keys() if fnmatch.fnmatch(p, item)]
                if matching_params:
                    return matching_params
                return []

        def resolve_group(entry, stage_name):
            params_list = []
            overrides = {}
            if isinstance(entry, str):
                tokens = entry.split()
                is_one_by_one = "one-by-one" in tokens
                base_tokens = [t for t in tokens if t != "one-by-one" and not t.startswith("wlmin=") and not t.startswith("wlmax=")]
                for t in tokens:
                    if t.startswith("wlmin="):
                        overrides['wlmin'] = float(t.split("=")[1])
                    elif t.startswith("wlmax="):
                        overrides['wlmax'] = float(t.split("=")[1])
                for item in base_tokens:
                    params_list.extend(resolve_single_param_or_group(item))
            elif isinstance(entry, list):
                is_one_by_one = "one-by-one" in entry
                for item in entry:
                    if item == "one-by-one" or isinstance(item, str) and (item.startswith("wlmin=") or item.startswith("wlmax=")):
                        if item.startswith("wlmin="):
                            overrides['wlmin'] = float(item.split("=")[1])
                        elif item.startswith("wlmax="):
                            overrides['wlmax'] = float(item.split("=")[1])
                        continue
                    params_list.extend(resolve_single_param_or_group(item))
            else:
                raise ValueError(f"Stage definition for '{stage_name}' must be a string or list")

            if is_one_by_one:
                sub_stages = []
                for i, param in enumerate(params_list):
                    var_part = param.split("_")[-1] if "_" in param else param
                    sub_name = f"{stage_name}_{var_part}" if len(params_list) > 1 else stage_name
                    sub_stages.append((sub_name, [param], overrides.copy()))
                return sub_stages
            return [(stage_name, params_list, overrides)]

        # Stages table with expanded stages
        stage_data = []
        for stage_name, stage_def in self.stages.items():
            resolved = resolve_group(stage_def, stage_name)
            for sub_name, params, overrides in resolved:
                param_str = ', '.join(params)
                if overrides:
                    param_str += f" (wlmin={overrides.get('wlmin', 'default')}, wlmax={overrides.get('wlmax', 'default')})"
                stage_data.append({
                    'Stage': sub_name,
                    'Parameters': param_str
                })
        stage_df = pd.DataFrame(stage_data)
        stage_html = stage_df.to_html(index=False, classes='table table-striped', border=0)

        html = f"""
        <div>
            <h4>TransmissionModel: {self.cross_section.name}</h4>
            <h5>Parameters</h5>
            {param_html}
            <h5>Fitting Stages</h5>
            {stage_html}
        </div>
        """
        return html
    def show_available_params(self, show_groups=True, show_params=True):
        """
        Display available parameter groups and individual parameters for Rietveld fitting.
        
        Parameters
        ----------
        show_groups : bool, optional
            If True, show predefined parameter groups
        show_params : bool, optional
            If True, show all individual parameters
        """
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
            "lattice": [p for p in self.params if p in ["a", "b", "c"]],
            "extinction": [p for p in self.params if p.startswith("ext_")],
            "sans": [p for p in self.params if p == "sans" or re.compile(r"sans\d+").match(p) or p.startswith("sans_")],
            "orientation": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ") or p.startswith("η")],
            "mosaicity": [p for p in self.params if p.startswith("η")],
            "thetas": [p for p in self.params if p.startswith("θ")],
            "phis": [p for p in self.params if p.startswith("ϕ")],
            "angles": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ")],
            "temperature": [p for p in ["temp"] if p in self.params],
        }
        if show_groups:
            print("Available parameter groups:")
            print("=" * 30)

            for group_name, params in group_map.items():
                if params:  # Only show groups with available parameters
                    print(f"  '{group_name}': {params}")
            
        if show_params:
            if show_groups:
                print("\nAll individual parameters:")
                print("=" * 30)
            else:
                print("Available parameters:")
                print("=" * 20)
                
            for param_name, param in self.params.items():
                vary_status = "vary" if param.vary else "fixed"
                print(f"  {param_name}: {param.value:.6g} ({vary_status})")
                
        print("\nExample usage:")
        print("=" * 15)
        print("# Using predefined groups:")
        print('param_groups = ["basic", "background", "extinction"]')
        print("\n# Using individual parameters:")
        print('param_groups = [["norm", "thickness"], ["b0", "ext_l2"]]')
        print("\n# Using named stages:")
        print('param_groups = {"scale": ["norm"], "sample": ["thickness", "extinction"]}')
        print("\n# Mixed approach:")
        print('param_groups = ["basic", ["b0", "ext_l2"], "lattice"]')
        print("\n# One-by-one expansion:")
        print('stages = {"angles_one": "angles one-by-one"}  # Expands to sub-stages for each angle')

    def set_cross_section(self, xs: 'CrossSection', inplace: bool = True) -> 'TransmissionModel':
        """
        Set a new cross-section for the model.

        Parameters
        ----------
        xs : CrossSection
            The new cross-section to apply.
        inplace : bool, optional
            If True, modify the current object. If False, return a new modified object, 
            by default True.

        Returns
        -------
        TransmissionModel
            The updated model (either modified in place or a new instance).
        """
        if inplace:
            self.cross_section = xs
            params = self._make_weight_params()
            self.params += params
            return self
        else:
            new_self = deepcopy(self)
            new_self.cross_section = xs
            params = new_self._make_weight_params()
            new_self.params += params
            return new_self
    def update_params(self, params: dict = {}, values_only: bool = True, inplace: bool = True):
        """
        Update the parameters of the model.

        Parameters
        ----------
        params : dict
            Dictionary of new parameters to update.
        values_only : bool, optional
            If True, update only the values of the parameters, by default True.
        inplace : bool, optional
            If True, modify the current object. If False, return a new modified object, 
            by default True.
        """
        if inplace:
            if values_only:
                for param in params:
                    self.params[param].set(value=params[param].value)
            else:
                self.params = params
        else:
            new_self = deepcopy(self)
            if values_only:
                for param in params:
                    new_self.params[param].set(value=params[param].value)
            else:
                new_self.params = params
            return new_self  # Ensure a return statement in the non-inplace scenario.
    def vary_all(self, vary: Optional[bool] = None, except_for: List[str] = []):
        """
        Toggle the 'vary' attribute for all model parameters.

        Parameters
        ----------
        vary : bool, optional
            The value to set for all parameters' 'vary' attribute.
        except_for : list of str, optional
            List of parameter names to exclude from this operation, by default [].
        """
        if vary is not None:
            for param in self.params:
                if param not in except_for:
                    self.params[param].set(vary=vary)
    def _tof_correction(self, E, L0: float = 1.0, t0: float = 0.0, **kwargs):
        """
        Apply a time-of-flight (TOF) correction to the energy values.

        Parameters
        ----------
        E : float or array-like
            The energy values to correct.
        L0 : float, optional
            The scale factor for the flight path, by default 1.0.
        t0 : float, optional
            The time offset for the correction, by default 0.0.
        kwargs : dict, optional
            Additional arguments (currently unused).

        Returns
        -------
        np.ndarray
            The corrected energy values.
        """
        tof = utils.energy2time(E, self.tof_length)
        dtof = (1.0 - L0) * tof + t0
        E = utils.time2energy(tof + dtof, self.tof_length)
        return E
    def group_weights(self, weights=None, vary=True, **groups):
        """
        Define softmax-normalized weight fractions for grouped phases, using shared `p1`, `p2`, ...
        parameters for internal group ratios, and global `group_<name>` parameters for relative group weights.

        Each group is normalized internally, and all groups sum to 1. Internal variation can be
        controlled per-group using the `vary` argument. Shared `pX` parameters are reused across groups.

        Parameters
        ----------
        weights : list of float, optional
            Initial relative weights between groups. Will be normalized. If not provided,
            all groups get equal initial weight.
        vary : bool or list of bool
            Whether to vary internal `pX` parameters of each group during fitting.
            Can be a single bool (applies to all groups), or a list of bools per group.
            Group-level weights always vary.
        **groups : dict[str, str | list[str]]
            Define each group by either:
            - a wildcard string (e.g., "inconel*")
            - or a list of phase names (e.g., ["inconel1", "inconel2"])

        Returns
        -------
        self : the model object

        Notes
        -----
        - This method reuses or creates global `p1`, `p2`, ... parameters to control phase weights.
        - Phase names are sanitized (dashes replaced with underscores).
        - The total sum of all phases will be 1.

        Examples
        --------
        >>> model = nbragg.TransmissionModel(xs)

        # Use wildcards and allow internal variation in both groups
        >>> model.group_weights(
        ...     inconel="inconel*",
        ...     steel="steel*",
        ...     weights=[0.7, 0.3],
        ...     vary=True
        ... )

        # Set internal variation only in 'inconel' group
        >>> model.group_weights(
        ...     inconel="inconel*",
        ...     steel="steel*",
        ...     weights=[0.5, 0.5],
        ...     vary=[True, False]
        ... )

        # Explicit group definitions (list of phases)
        >>> model.group_weights(
        ...     powder=["inconel0", "inconel1", "steel_powder"],
        ...     bulk=["steel0", "steel1", "steel2"],
        ...     weights=[0.2, 0.8],
        ...     vary=False
        ... )
        """
        import fnmatch
        from numpy import log
        import lmfit

        self.params = getattr(self, "params", lmfit.Parameters())
        all_phases = list(self._materials.keys())
        group_names = list(groups.keys())
        num_groups = len(group_names)

        # Normalize 'vary'
        if isinstance(vary, bool):
            vary = [vary] * num_groups
        assert len(vary) == num_groups, "Length of `vary` must match number of groups"

        # Normalize 'weights'
        if weights is None:
            weights = [1.0] * num_groups
        assert len(weights) == num_groups, "Length of `weights` must match number of groups"

        # Resolve wildcard groups
        resolved_groups = {}
        for name, spec in groups.items():
            if isinstance(spec, str):
                matched = sorted(fnmatch.filter(all_phases, spec))
            elif isinstance(spec, list):
                matched = spec
            else:
                raise ValueError(f"Group '{name}' must be a string or list of phase names")
            if not matched:
                raise ValueError(f"No phases matched for group '{name}' using '{spec}'")
            resolved_groups[name] = matched

        # Add group weight softmax parameters: g1, g2, ...
        for i in range(num_groups - 1):
            val = log(weights[i] / weights[-1])
            self.params.add(f"g{i+1}", value=val, min=-14, max=14, vary=True)

        denom = " + ".join([f"exp(g{i+1})" for i in range(num_groups - 1)] + ["1"])
        for i, group in enumerate(group_names[:-1]):
            self.params.add(f"group_{group}", expr=f"exp(g{i+1}) / ({denom})")
        self.params.add(f"group_{group_names[-1]}", expr=f"1 / ({denom})")

        # Clear any existing p-parameters that might conflict
        # We'll rebuild them from scratch
        existing_p_params = [name for name in self.params.keys() if name.startswith('p') and name[1:].isdigit()]
        for p_name in existing_p_params:
            del self.params[p_name]
        
        # Clear any existing phase parameters that will be rebuilt
        all_group_phases = []
        for phases in resolved_groups.values():
            all_group_phases.extend([phase.replace("-", "") for phase in phases])
        
        for phase_name in all_group_phases:
            if phase_name in self.params:
                del self.params[phase_name]

        # Assign p1, p2, ..., shared across all groups — exactly N-1 per group
        p_index = 1

        for group_i, group_name in enumerate(group_names):
            phases = resolved_groups[group_name]
            group_frac = f"group_{group_name}"
            N = len(phases)

            if N == 1:
                phase_clean = phases[0].replace("-", "")
                self.params.add(phase_clean, expr=group_frac)
                continue

            # Create exactly N-1 parameters for this group
            group_pnames = []
            for i in range(N - 1):  # Only N−1 softmax params per group
                pname = f"p{p_index}"
                p_index += 1
                
                # Get initial value from material weights
                phase = phases[i]
                val = log(self._materials[phase]["weight"] / self._materials[phases[-1]]["weight"])
                
                # Add the parameter if it doesn't exist, or update vary if it does
                if pname in self.params:
                    self.params[pname].set(vary=vary[group_i])
                else:
                    self.params.add(pname, value=val, min=-14, max=14, vary=vary[group_i])
                
                group_pnames.append(pname)

            # Build denominator expression
            denom_terms = [f"exp({pname})" for pname in group_pnames]
            denom_expr = "1 + " + " + ".join(denom_terms)

            # Add expressions for first N-1 phases
            for i, phase in enumerate(phases[:-1]):
                phase_clean = phase.replace("-", "")
                pname = group_pnames[i]
                self.params.add(phase_clean, expr=f"{group_frac} * exp({pname}) / ({denom_expr})")

            # Add expression for the last phase (reference phase)
            final_phase = phases[-1].replace("-", "")
            self.params.add(final_phase, expr=f"{group_frac} / ({denom_expr})")
