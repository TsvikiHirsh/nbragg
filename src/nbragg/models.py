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
import ipywidgets as widgets
from IPython.display import display
from matplotlib.patches import Rectangle
import fnmatch
import re
from numpy import log


class TransmissionModel(lmfit.Model):
    def __init__(self, cross_section,
                params: "lmfit.Parameters" = None,
                response: str = "jorgensen",
                background: str = "polynomial3",
                tof_length: float = 9,
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
        cross_section : callable
            A function that takes energy (E) as input and returns the cross section.
        response : str, optional
            The type of response function to use, by default "jorgensen".
        background : str, optional
            The type of background function to use, by default "polynomial3".
        tof_length : float, optional
            The flight path length in [m]
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
        """
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
            "basic": True,  # Always include basic parameters
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

        .. math:: T(\lambda) = \text{norm} \cdot e^{- \sigma \cdot \text{thickness} \cdot n} \cdot (1 - \text{bg}) + \text{bg}
        
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
    
    def fit(self, data, params=None, wlmin: float = 1., wlmax: float = 6.,
            method: str = "rietveld",
            xtol: float = None, ftol: float = None, gtol: float = None,
            verbose: bool = False,
            progress_bar: bool = True,
            stages: Optional[Union[str, Dict[str, Union[str, List[str]]]]] = None,
            **kwargs):
        """
        Fit the model to data.

        This method supports multiple fitting approaches:
        - **Standard single-stage fitting** (`method="least-squares"`)
        - **True Rietveld-style refinement** (`method="rietveld"`, default) - parameters accumulate across stages
        - **Staged sequential refinement** (`method="staged"`) - parameters are frozen after each stage

        Parameters
        ----------
        data : pandas.DataFrame or Data or array-like
            The input data.
        params : lmfit.Parameters, optional
            Parameters to use for fitting. If None, uses the model's default parameters.
        wlmin, wlmax : float, optional
            Minimum and maximum wavelength for fitting.
        method : str, optional
            Fitting method: "least-squares", "rietveld", or "staged" (default is "rietveld").
        xtol, ftol, gtol : float, optional
            Convergence tolerances (passed to `lmfit`).
        verbose : bool, optional
            If True, prints detailed fitting information.
        progress_bar : bool, optional
            If True, shows a progress bar for fitting.
        stages : str or dict, optional
            Fitting stages. Can be "all" or a dictionary of stage definitions.
            If None, uses self.stages.
        **kwargs
            Additional keyword arguments passed to `lmfit.Model.fit`.

        Returns
        -------
        lmfit.model.ModelResult
            The fit result object.

        Examples
        --------
        >>> import nbragg
        >>> # Create a sample cross-section, data and model
        >>> xs = nbragg.CrossSection(...)  # Assume a valid CrossSection
        >>> data = nbragg.Data(...)  # Assume valid Data
        >>> model = nbragg.TransmissionModel(xs, vary_background=True, vary_weights=True)

        # Default Rietveld fitting with automatic stages
        >>> result = model.fit(data)
        
        # Single-stage fitting with all vary=True parameters
        >>> result = model.fit(data, stages="all")
        
        # Custom stages for Rietveld fitting
        >>> stages = {"background": "background", "scale": ["norm", "thickness"]}
        >>> result = model.fit(data, stages=stages)
        
        # Set custom stages on the model and fit
        >>> model.stages = {"stage1": ["b0", "b1"], "stage2": "all"}
        >>> result = model.fit(data)
        """
        # Handle stages argument
        if stages is not None:
            if isinstance(stages, str) and stages == "all":
                stages = {"all": "all"}
            elif not isinstance(stages, dict):
                raise ValueError("Stages must be 'all' or a dictionary")
        else:
            stages = self.stages

        # Route to multi-stage fitting if requested
        if method in ["rietveld", "staged"]:
            return self._multistage_fit(
                data, params, wlmin, wlmax,
                method=method,
                verbose=verbose,
                progress_bar=progress_bar,
                stages=stages,
                **kwargs
            )

        # Prepare fit kwargs
        fit_kws = kwargs.pop("fit_kws", {})
        if xtol is not None: fit_kws.setdefault("xtol", xtol)
        if ftol is not None: fit_kws.setdefault("ftol", ftol)
        if gtol is not None: fit_kws.setdefault("gtol", gtol)
        kwargs["fit_kws"] = fit_kws

        # Try tqdm for progress
        try:
            from tqdm.notebook import tqdm
        except ImportError:
            from tqdm.auto import tqdm
            

        # If progress_bar=True, wrap the fit in tqdm
        if progress_bar:
            pbar = tqdm(total=1, desc="Fitting", disable=not progress_bar)
        else:
            pbar = None

        # Prepare input data
        if isinstance(data, pandas.DataFrame):
            data = data.query(f"{wlmin} < wavelength < {wlmax}")
            weights = kwargs.get("weights", 1. / data["err"].values)
            fit_result = super().fit(
                data["trans"].values,
                params=params or self.params,
                weights=weights,
                wl=data["wavelength"].values,
                method=method,
                **kwargs
            )

        elif isinstance(data, Data):
            data = data.table.query(f"{wlmin} < wavelength < {wlmax}")
            weights = kwargs.get("weights", 1. / data["err"].values)
            fit_result = super().fit(
                data["trans"].values,
                params=params or self.params,
                weights=weights,
                wl=data["wavelength"].values,
                method=method,
                **kwargs
            )

        else:
            fit_result = super().fit(
                data,
                params=params or self.params,
                method=method,
                **kwargs
            )

        if pbar:
            pbar.set_postfix({"redchi": f"{fit_result.redchi:.4g}"})
            pbar.update(1)
            pbar.close()

        # Attach results
        self.fit_result = fit_result
        fit_result.plot = self.plot
        fit_result.plot_total_xs = self.plot_total_xs
        fit_result.show_available_params = self.show_available_params

        if self.response is not None:
            fit_result.response = self.response
            fit_result.response.params = fit_result.params
        if self.background is not None:
            fit_result.background = self.background

        return fit_result
    
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


    def _get_stage_parameters(self, stage_def: Union[str, List[str]]) -> List[str]:
        """Helper method to get parameters associated with a stage definition."""
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
        if stage_def == "all":
            return [p for p in self.params if self.params[p].vary]
        if isinstance(stage_def, str):
            return group_map.get(stage_def, [stage_def] if stage_def in self.params else [])
        params = []
        for item in stage_def:
            if item in group_map:
                params.extend(group_map[item])
            elif item in self.params:
                params.append(item)
            else:
                matching_params = [p for p in self.params.keys() if fnmatch.fnmatch(p, item)]
                params.extend(matching_params)
        return list(dict.fromkeys(params))  # Remove duplicates while preserving order

    def _multistage_fit(self, data, params: "lmfit.Parameters" = None, wlmin: float = 1, wlmax: float = 8,
                        method: str = "staged",
                        verbose=False, progress_bar=True,
                        stages=None,
                        **kwargs):
        """ 
        Perform multi-stage fitting with two different strategies:
        
        - "rietveld": True Rietveld refinement where parameters accumulate across stages
        - "staged": Sequential staged refinement where parameters are frozen after each stage
        
        Parameters
        ----------
        data : pandas.DataFrame or Data
            The input data containing wavelength and transmission values.
        params : lmfit.Parameters, optional
            Initial parameters for the fit. If None, uses the model's default parameters.
        wlmin : float, optional default=1
            Default minimum wavelength for fitting.
        wlmax : float, optional default=8
            Default maximum wavelength for fitting.
        method : str, optional
            Fitting method: "rietveld" or "staged".
        verbose : bool, optional
            If True, prints detailed information about each fitting stage.
        progress_bar : bool, optional
            If True, shows a progress bar for each fitting stage.
        stages : dict, optional
            Dictionary of stage definitions. If None, uses self.stages.
        **kwargs
            Additional keyword arguments for the fit method.

        Returns
        -------
        fit_result : lmfit.ModelResult
            The final fit result after all stages.
        """
        from copy import deepcopy
        import sys
        import warnings
        import re
        import fnmatch
        import pandas
        try:
            from tqdm.notebook import tqdm
        except ImportError:
            from tqdm.auto import tqdm
        import pickle

        if method not in ["rietveld", "staged"]:
            raise ValueError(f"Invalid multi-stage method: {method}. Use 'rietveld' or 'staged'.")

        # User-friendly group name mapping
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
            """Resolve a single parameter name or group name to a list of parameters."""
            if item == "all":
                return [p for p in self.params if self.params[p].vary]
            elif item == "one-by-one":
                return []  # Handled separately in resolve_group
            elif item in group_map:
                resolved = group_map[item]
                if verbose:
                    print(f"  Resolved group '{item}' to: {resolved}")
                return resolved
            elif item in self.params:
                if verbose:
                    print(f"  Found parameter: {item}")
                return [item]
            else:
                matching_params = [p for p in self.params.keys() if fnmatch.fnmatch(p, item)]
                if matching_params:
                    if verbose:
                        print(f"  Pattern '{item}' matched: {matching_params}")
                    return matching_params
                else:
                    warnings.warn(f"Unknown parameter or group: '{item}'. Available parameters: {list(self.params.keys())}")
                    return []

        def resolve_group(entry, stage_name):
            """
            Resolve a group entry to a list of parameters and overrides.
            If "one-by-one" is detected in the entry list, expand all parameters into sub-stages.
            """
            if isinstance(entry, str):
                tokens = entry.split()
                params_list = []
                overrides = {}
                is_one_by_one = "one-by-one" in tokens
                if is_one_by_one:
                    idx = tokens.index("one-by-one")
                    base_tokens = tokens[:idx]
                    post_tokens = tokens[idx + 1:]
                    base_entry = " ".join(base_tokens)
                    # Process base for params
                    base_items = base_entry.split() if base_entry else []
                    for it in base_items:
                        params_list.extend(resolve_single_param_or_group(it))
                    # Process post for overrides
                    for tok in post_tokens:
                        if tok.startswith("wlmin="):
                            k, v = tok.split("=")
                            overrides['wlmin'] = float(v)
                        elif tok.startswith("wlmax="):
                            k, v = tok.split("=")
                            overrides['wlmax'] = float(v)
                else:
                    # Normal processing
                    for it in tokens:
                        if it.startswith("wlmin="):
                            k, v = it.split("=")
                            overrides['wlmin'] = float(v)
                        elif it.startswith("wlmax="):
                            k, v = it.split("=")
                            overrides['wlmax'] = float(v)
                        else:
                            params_list.extend(resolve_single_param_or_group(it))
                if is_one_by_one:
                    sub_stages = []
                    for i, param in enumerate(params_list):
                        var_part = param.split("_")[-1] if "_" in param else param
                        sub_name = f"{stage_name}_{var_part}" if len(params_list) > 1 else stage_name
                        sub_stages.append((sub_name, [param], overrides.copy()))
                    return sub_stages
                return [(stage_name, params_list, overrides)]
            elif isinstance(entry, list):
                params_list = []
                overrides = {}
                is_one_by_one = "one-by-one" in entry
                for item in entry:
                    if item == "one-by-one":
                        continue
                    if isinstance(item, str) and item.startswith("wlmin="):
                        try:
                            overrides['wlmin'] = float(item.split("=", 1)[1])
                            if verbose:
                                print(f"  Override wlmin detected: {overrides['wlmin']}")
                        except ValueError:
                            warnings.warn(f"Invalid wlmin value in group: {item}")
                    elif isinstance(item, str) and item.startswith("wlmax="):
                        try:
                            overrides['wlmax'] = float(item.split("=", 1)[1])
                            if verbose:
                                print(f"  Override wlmax detected: {overrides['wlmax']}")
                        except ValueError:
                            warnings.warn(f"Invalid wlmax value in group: {item}")
                    else:
                        params_list.extend(resolve_single_param_or_group(item))
                if is_one_by_one:
                    sub_stages = []
                    for i, param in enumerate(params_list):
                        var_part = param.split("_")[-1] if "_" in param else param
                        sub_name = f"{stage_name}_{var_part}" if len(params_list) > 1 else stage_name
                        sub_stages.append((sub_name, [param], overrides.copy()))
                    return sub_stages
                return [(stage_name, params_list, overrides)]
            else:
                raise ValueError(f"Stage definition for '{stage_name}' must be a string or list, got {type(entry)}")

        # Handle stages input
        expanded_stages = []
        if isinstance(stages, dict):
            for stage_name, entry in stages.items():
                resolved = resolve_group(entry, stage_name)
                expanded_stages.extend(resolved)
        else:
            raise ValueError("Stages must be a dictionary")

        # Remove any empty stages
        filtered = [(n, g, o) for n, g, o in zip(*zip(*expanded_stages)) if g]
        if not filtered:
            raise ValueError("No valid stages found. Check your stage definitions.")
        stage_names, resolved_stages, stage_overrides = zip(*filtered)

        if verbose:
            refinement_type = "True Rietveld (accumulative)" if method == "rietveld" else "Staged sequential"
            print(f"\n{refinement_type} fitting stages with possible wavelength overrides:")
            for i, (name, group, ov) in enumerate(zip(stage_names, resolved_stages, stage_overrides)):
                print(f"  {name}: {group if group else 'all vary=True parameters'}  overrides: {ov}")

        # Store for summary or introspection
        self._stage_param_groups = list(resolved_stages)
        self._stage_names = list(stage_names)
        self._fitting_method = method

        params = deepcopy(params or self.params)

        # Setup tqdm iterator
        try:
            from tqdm.notebook import tqdm
            if 'ipykernel' in sys.modules:
                iterator = tqdm(
                    zip(stage_names, resolved_stages, stage_overrides),
                    desc=f"{'Rietveld' if method == 'rietveld' else 'Staged'} Fit",
                    disable=not progress_bar,
                    total=len(stage_names)
                )
            else:
                iterator = tqdm(
                    zip(stage_names, resolved_stages, stage_overrides),
                    desc=f"{'Rietveld' if method == 'rietveld' else 'Staged'} Fit",
                    disable=not progress_bar,
                    total=len(stage_names)
                )
        except ImportError:
            iterator = tqdm(
                zip(stage_names, resolved_stages, stage_overrides),
                desc=f"{'Rietveld' if method == 'rietveld' else 'Staged'} Fit",
                disable=not progress_bar,
                total=len(stage_names)
            )

        stage_results = []
        stage_summaries = []
        cumulative_params = set()  # Track parameters that have been refined (for rietveld method)

        def extract_pickleable_attributes(fit_result):
            safe_attrs = [
                'params', 'success', 'residual', 'chisqr', 'redchi', 'aic', 'bic',
                'nvarys', 'ndata', 'nfev', 'message', 'lmdif_message', 'cov_x',
                'method', 'flatchain', 'errorbars', 'ci_out'
            ]

            class PickleableResult:
                pass

            result = PickleableResult()

            for attr in safe_attrs:
                if hasattr(fit_result, attr):
                    try:
                        value = getattr(fit_result, attr)
                        pickle.dumps(value)
                        setattr(result, attr, value)
                    except (TypeError, ValueError, AttributeError):
                        if verbose:
                            print(f"Skipping non-pickleable attribute: {attr}")
                        continue

            return result

        for stage_idx, (stage_name, group, overrides) in enumerate(iterator):
            stage_num = stage_idx + 1

            # Use overrides or fallback to global wlmin, wlmax
            stage_wlmin = overrides.get('wlmin', wlmin)
            stage_wlmax = overrides.get('wlmax', wlmax)

            if verbose:
                group_display = group if group else "all vary=True parameters"
                print(f"\n{stage_name}: Fitting parameters {group_display} with wavelength range [{stage_wlmin}, {stage_wlmax}]")

            # Filter data for this stage
            if isinstance(data, pandas.DataFrame):
                stage_data = data.query(f"{stage_wlmin} < wavelength < {stage_wlmax}")
                wavelengths = stage_data["wavelength"].values
                trans = stage_data["trans"].values
                weights = kwargs.get("weights", 1. / stage_data["err"].values)
            elif isinstance(data, Data):
                stage_data = data.table.query(f"{stage_wlmin} < wavelength < {stage_wlmax}")
                wavelengths = stage_data["wavelength"].values
                trans = stage_data["trans"].values
                weights = kwargs.get("weights", 1. / stage_data["err"].values)
            else:
                raise ValueError("Multi-stage fitting requires wavelength-based input data.")

            # Set parameter vary status based on method
            if method == "rietveld":
                # True Rietveld: accumulate parameters across stages
                cumulative_params.update(group if group else [p for p in self.params if self.params[p].vary])
                
                # Freeze all parameters first
                for p in params.values():
                    p.vary = False
                
                # Unfreeze all parameters that have been introduced so far
                unfrozen_count = 0
                for name in cumulative_params:
                    if name in params:
                        params[name].vary = True
                        unfrozen_count += 1
                        if verbose and (name in group or not group):
                            print(f"  New parameter: {name}")
                        elif verbose:
                            print(f"  Continuing: {name}")
                    else:
                        if name in group or not group:  # Only warn for new parameters
                            warnings.warn(f"Parameter '{name}' not found in params")
                
                if verbose:
                    print(f"  Total active parameters: {unfrozen_count}")
                    
            elif method == "staged":
                # Staged: only current group parameters vary
                for p in params.values():
                    p.vary = False

                unfrozen_count = 0
                active_params = group if group else [p for p in self.params if self.params[p].vary]
                for name in active_params:
                    if name in params:
                        params[name].vary = True
                        unfrozen_count += 1
                        if verbose:
                            print(f"  Unfrozen: {name}")
                    else:
                        warnings.warn(f"Parameter '{name}' not found in params")

            if unfrozen_count == 0:
                warnings.warn(f"No parameters were unfrozen in {stage_name}. Skipping this stage.")
                continue

            # Perform fitting
            try:
                fit_result = super().fit(
                    trans,
                    params=params,
                    wl=wavelengths,
                    weights=weights,
                    method="leastsq",
                    **kwargs
                )
            except Exception as e:
                warnings.warn(f"Fitting failed in {stage_name}: {e}")
                continue

            # Extract pickleable part
            stripped_result = extract_pickleable_attributes(fit_result)

            stage_results.append(stripped_result)

            # Build summary
            if method == "rietveld":
                varied_params = list(cumulative_params)
            else:
                varied_params = group if group else [p for p in self.params if self.params[p].vary]
                
            summary = {
                "stage": stage_num,
                "stage_name": stage_name,
                "fitted_params": group if group else ["all vary=True"],
                "active_params": varied_params,
                "wlmin": stage_wlmin,
                "wlmax": stage_wlmax,
                "redchi": fit_result.redchi,
                "method": method
            }
            for name, par in fit_result.params.items():
                summary[f"{name}_value"] = par.value
                summary[f"{name}_stderr"] = par.stderr
                summary[f"{name}_vary"] = name in varied_params
            stage_summaries.append(summary)

            method_display = "Rietveld" if method == "rietveld" else "Staged"
            iterator.set_description(f"{method_display} {stage_num}/{len(stage_names)}")
            iterator.set_postfix({"stage": stage_name, "reduced χ²": f"{fit_result.redchi:.4g}"})

            # Update params for next stage
            params = fit_result.params

            if verbose:
                print(f"  {stage_name} completed. χ²/dof = {fit_result.redchi:.4f}")

        if not stage_results:
            raise RuntimeError("No successful fitting stages completed")

        self.fit_result = fit_result
        self.fit_stages = stage_results
        self.stages_summary = self._create_stages_summary_table_enhanced(
            stage_results, resolved_stages, stage_names, method=method
        )

        # Attach plotting methods and other attributes
        fit_result.plot = self.plot
        fit_result.plot_total_xs = self.plot_total_xs
        fit_result.plot_stage_progression = self.plot_stage_progression
        fit_result.plot_chi2_progression = self.plot_chi2_progression
        if self.response is not None:
            fit_result.response = self.response
            fit_result.response.params = fit_result.params
        if self.background is not None:
            fit_result.background = self.background

        fit_result.stages_summary = self.stages_summary
        fit_result.show_available_params = self.show_available_params
        return fit_result




    def _create_stages_summary_table_enhanced(self, stage_results, resolved_param_groups, stage_names=None, 
                                            method="rietveld", color=True):
        import pandas as pd
        import numpy as np

        # --- Build the DataFrame ---
        all_param_names = list(stage_results[-1].params.keys())
        stage_data = {}
        if stage_names is None:
            stage_names = [f"Stage_{i+1}" for i in range(len(stage_results))]

        cumulative_params = set()  # Track cumulative parameters for Rietveld method

        for stage_idx, stage_result in enumerate(stage_results):
            stage_col = stage_names[stage_idx] if stage_idx < len(stage_names) else f"Stage_{stage_idx + 1}"
            stage_data[stage_col] = {'value': {}, 'stderr': {}, 'vary': {}}
            
            # Determine which parameters varied in this stage
            if method == "rietveld":
                # For Rietveld: accumulate parameters
                cumulative_params.update(resolved_param_groups[stage_idx])
                varied_in_stage = cumulative_params.copy()
            else:
                # For staged: only current group
                varied_in_stage = set(resolved_param_groups[stage_idx])

            for param_name in all_param_names:
                if param_name in stage_result.params:
                    param = stage_result.params[param_name]
                    stage_data[stage_col]['value'][param_name] = param.value
                    stage_data[stage_col]['stderr'][param_name] = param.stderr if param.stderr is not None else np.nan
                    stage_data[stage_col]['vary'][param_name] = param_name in varied_in_stage
                else:
                    stage_data[stage_col]['value'][param_name] = np.nan
                    stage_data[stage_col]['stderr'][param_name] = np.nan
                    stage_data[stage_col]['vary'][param_name] = False

            redchi = stage_result.redchi if hasattr(stage_result, 'redchi') else np.nan
            stage_data[stage_col]['value']['redchi'] = redchi
            stage_data[stage_col]['stderr']['redchi'] = np.nan
            stage_data[stage_col]['vary']['redchi'] = np.nan

        # Create DataFrame
        data_for_df = {}
        for stage_col in stage_data:
            for metric in ['value', 'stderr', 'vary']:
                data_for_df[(stage_col, metric)] = stage_data[stage_col][metric]

        df = pd.DataFrame(data_for_df)
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=['Stage', 'Metric'])
        all_param_names_with_redchi = all_param_names + ['redchi']
        df = df.reindex(all_param_names_with_redchi)

        # --- Add initial values column ---
        initial_values = {}
        for param_name in all_param_names:
            initial_values[param_name] = self.params[param_name].value if param_name in self.params else np.nan
        initial_values['redchi'] = np.nan

        initial_df = pd.DataFrame({('Initial', 'value'): initial_values})
        df = pd.concat([initial_df, df], axis=1)

        if not color:
            return df

        styler = df.style

        # 1) Highlight vary=True cells with different colors for different methods
        vary_cols = [col for col in df.columns if col[1] == 'vary']
        if method == "rietveld":
            # Light green for Rietveld (accumulative)
            def highlight_vary_rietveld(s):
                return ['background-color: lightgreen' if v is True else '' for v in s]
            for col in vary_cols:
                styler = styler.apply(highlight_vary_rietveld, subset=[col], axis=0)
        else:
            # Light blue for staged (sequential)
            def highlight_vary_staged(s):
                return ['background-color: lightblue' if v is True else '' for v in s]
            for col in vary_cols:
                styler = styler.apply(highlight_vary_staged, subset=[col], axis=0)

        # 2) Highlight redchi row's value cells (moccasin)
        def highlight_redchi_row(row):
            if row.name == 'redchi':
                return ['background-color: moccasin' if col[1] == 'value' else '' for col in df.columns]
            return ['' for _ in df.columns]
        styler = styler.apply(highlight_redchi_row, axis=1)

        # 3) Highlight value cells by fractional change with red hues (ignore <1%)
        value_cols = [col for col in df.columns if col[1] == 'value']

        # Calculate % absolute change between consecutive columns (Initial → Stage1 → Stage2 ...)
        changes = pd.DataFrame(index=df.index, columns=value_cols, dtype=float)
        prev_col = None
        for col in value_cols:
            if prev_col is None:
                # No previous for initial column, so zero changes here
                changes[col] = 0.0
            else:
                prev_vals = df[prev_col].astype(float)
                curr_vals = df[col].astype(float)
                with np.errstate(divide='ignore', invalid='ignore'):
                    pct_change = np.abs((curr_vals - prev_vals) / prev_vals) * 100
                pct_change = pct_change.replace([np.inf, -np.inf], np.nan).fillna(0.0)
                changes[col] = pct_change
            prev_col = col

        max_change = changes.max().max()
        # Normalize by max change, to get values in [0,1]
        norm_changes = changes / max_change if max_change > 0 else changes

        def red_color(val):
            # Ignore changes less than 1%
            if pd.isna(val) or val < 1:
                return ''
            # val in [0,1], map to red intensity
            # 0 -> white (255,255,255)
            # 1 -> dark red (255,100,100)
            r = 255
            g = int(255 - 155 * val)
            b = int(255 - 155 * val)
            return f'background-color: rgb({r},{g},{b})'

        for col in value_cols:
            styler = styler.apply(lambda s: [red_color(v) for v in norm_changes[col]], subset=[col], axis=0)

        return styler





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

    def plot(self, data=None, plot_bg: bool = True,    
            plot_dspace: bool = False, dspace_min: float = 1,    
            dspace_label_pos: float = 0.99, stage: int = None, **kwargs):    
        """    
        Plot the results of the fit or model.    
            
        Parameters    
        ----------    
        data : object, optional    
            Data object to show alongside the model (useful before performing the fit).    
            Should have wavelength, transmission, and error data accessible.    
        plot_bg : bool, optional    
            Whether to include the background in the plot, by default True.    
        plot_dspace: bool, optional    
            If True plots the 2*dspace and labels of that material that are larger than dspace_min    
        dspace_min: float, optional    
            The minimal dspace from which to plot the dspacing*2 lines    
        dspace_label_pos: float, optional    
            The position on the y-axis to plot the dspace label, e.g. 1 is at the top of the figure    
        stage: int, optional    
            If provided, plot results from a specific Rietveld fitting stage (1-indexed).    
            Only works if Rietveld fitting has been performed.    
        kwargs : dict, optional    
            Additional plot settings like color, marker size, etc.    
                
        Returns    
        -------    
        matplotlib.axes.Axes    
            The axes of the plot.    
                
        Notes    
        -----    
        This function generates a plot showing the transmission data, the best-fit curve,    
        and residuals. If `plot_bg` is True, it will also plot the background function.    
        Can be used both after fitting (using fit_result) or before fitting (using model params).    
        """    
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig, ax = plt.subplots(2, 1, sharex=True, height_ratios=[3.5, 1], figsize=(6, 5))    
            
        # Determine which results to use
        if stage is not None and hasattr(self, "fit_stages") and self.fit_stages:
            # Use specific stage results
            if stage < 1 or stage > len(self.fit_stages):
                raise ValueError(f"Stage {stage} not available. Available stages: 1-{len(self.fit_stages)}")
            
            # Get stage results
            stage_result = self.fit_stages[stage - 1]  # Convert to 0-indexed
            
            # We need to reconstruct the fit data from the original fit
            if hasattr(self, "fit_result") and self.fit_result is not None:
                wavelength = self.fit_result.userkws["wl"]    
                data_values = self.fit_result.data    
                err = 1. / self.fit_result.weights    
            else:
                raise ValueError("Cannot plot stage results without original fit data")
                
            # Use stage parameters to evaluate model
            params = stage_result.params
            best_fit = self.eval(params=params, wl=wavelength)
            residual = (data_values - best_fit) / err
            chi2 = stage_result.redchi if hasattr(stage_result, 'redchi') else np.sum(residual**2) / (len(data_values) - len(params))
            fit_label = f"Stage {stage} fit"
            
        elif hasattr(self, "fit_result") and self.fit_result is not None:    
            # Use final fit results    
            wavelength = self.fit_result.userkws["wl"]    
            data_values = self.fit_result.data    
            err = 1. / self.fit_result.weights    
            best_fit = self.fit_result.best_fit    
            residual = self.fit_result.residual    
            params = self.fit_result.params    
            chi2 = self.fit_result.redchi    
            fit_label = "Best fit"    
        else:    
            # Use model (no fit yet)    
            fit_label = "Model"    
            params = self.params  # Assuming model has params attribute    
                
            if data is not None:    
                # Extract data from provided data object    
                wavelength = data.table.wavelength    
                data_values = data.table.trans    
                err = data.table.err    
                    
                # Evaluate model at data wavelengths    
                best_fit = self.eval(params=params, wl=wavelength)    
                residual = (data_values - best_fit) / err    
                    
                # Calculate chi2 for the model    
                chi2 = np.sum(((data_values - best_fit) / err) ** 2) / (len(data_values) - len(params))    
            else:    
                # No data provided, just show model over some wavelength range    
                wavelength = np.linspace(1.0, 10.0, 1000)  # Adjust range as needed    
                data_values = np.nan * np.ones_like(wavelength)    
                err = np.nan * np.ones_like(wavelength)    
                best_fit = self.eval(params=params, wl=wavelength)    
                residual = np.nan * np.ones_like(wavelength)    
                chi2 = np.nan    
            
        # Plot settings    
        color = kwargs.pop("color", "seagreen")    
        ecolor = kwargs.pop("ecolor", "0.8")    
        title = kwargs.pop("title", self.cross_section.name)    
        ms = kwargs.pop("ms", 2)    
            
        # Plot data and best-fit/model    
        ax[0].errorbar(wavelength, data_values, err, marker="o", color=color, ms=ms,     
                    zorder=-1, ecolor=ecolor, label="Data")    
        ax[0].plot(wavelength, best_fit, color="0.2", label=fit_label)    
        ax[0].set_ylabel("Transmission")    
        ax[0].set_title(title)    
            
        # Plot residuals    
        ax[1].plot(wavelength, residual, color=color)    
        ax[1].set_ylabel("Residuals [1σ]")    
        ax[1].set_xlabel("λ [Å]")    
            
        # Plot background if requested    
        if plot_bg and self.background:    
            self.background.plot(wl=wavelength, ax=ax[0], params=params, **kwargs)    
            legend_labels = [fit_label, "Background", "Data"]    
        else:    
            legend_labels = [fit_label, "Data"]    
            
        # Set legend with chi2 value    
        ax[0].legend(legend_labels, fontsize=9, reverse=True,     
                    title=f"χ$^2$: {chi2:.2f}" if not np.isnan(chi2) else "χ$^2$: N/A")    
            
        # Plot d-spacing lines if requested    
        if plot_dspace:    
            for phase in self.cross_section.phases_data:    
                try:    
                    hkls = self.cross_section.phases_data[phase].info.hklList()    
                except:    
                    continue    
                for hkl in hkls:    
                    hkl = hkl[:3]    
                    dspace = self.cross_section.phases_data[phase].info.dspacingFromHKL(*hkl)    
                    if dspace >= dspace_min:    
                        trans = ax[0].get_xaxis_transform()    
                        ax[0].axvline(dspace*2, lw=1, color="0.4", zorder=-1, ls=":")    
                        if len(self.cross_section.phases) > 1:    
                            ax[0].text(dspace*2, dspace_label_pos, f"{phase} {hkl}",     
                                    color="0.2", zorder=-1, fontsize=8, transform=trans,     
                                    rotation=90, va="top", ha="right")    
                        else:    
                            ax[0].text(dspace*2, dspace_label_pos, f"{hkl}",     
                                    color="0.2", zorder=-1, fontsize=8, transform=trans,     
                                    rotation=90, va="top", ha="right")    
            
        plt.subplots_adjust(hspace=0.05)    
        return ax    

    def _make_basic_params(self):
        params = lmfit.Parameters()
        params.add("thickness", value=1., min=0.)
        params.add("norm", value=1., min=0.)
        return params

    def _make_temperature_params(self):
        params = lmfit.Parameters()
        params.add("temp", value=293.15, min=0.)  # Default temperature in Kelvin
        return params

    def _make_weight_params(self, vary=False):
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
        params = lmfit.Parameters()
        for phase in self.cross_section.phases:
            params.add(f"θ_{phase}", value=0., vary=vary)
            params.add(f"ϕ_{phase}", value=0., vary=vary)
            params.add(f"η_{phase}", value=0., min=0., vary=vary)
        return params

    def _tof_correction(self, E, **kwargs):
        L0 = kwargs.get("L0", self.tof_length)
        t0 = kwargs.get("t0", 0.)
        # Assuming energy correction based on TOF
        return E * (L0 / self.tof_length) + t0

    def _make_tof_params(self, vary=False, **kwargs):
        params = lmfit.Parameters()
        params.add("L0", value=1., min=0., max = 2., vary=vary)
        params.add("t0", value=0., vary=vary)
        return params

    def plot_total_xs(self, plot_bg: bool = True,     
                    plot_dspace: bool = False,     
                    dspace_min: float = 1,     
                    dspace_label_pos: float = 0.99,     
                    stage: int = None,
                    **kwargs):    
        """    
        Plot the results of the total cross-section fit.    

        Parameters    
        ----------    
        plot_bg : bool, optional    
            Whether to include the background in the plot, by default True.    
        plot_dspace: bool, optional    
            If True plots the 2*dspace and labels of that material that are larger than dspace_min    
        dspace_min: float, optional    
            The minimal dspace from which to plot the dspacing*2 lines    
        dspace_label_pos: float, optional    
            The position on the y-axis to plot the dspace label, e.g. 1 is at the top of the figure    
        stage: int, optional    
            If provided, plot results from a specific Rietveld fitting stage (1-indexed).    
            Only works if Rietveld fitting has been performed.    
        kwargs : dict, optional    
            Additional plot settings like color, marker size, etc.    
                
        Returns    
        -------    
        matplotlib.axes.Axes    
            The axes of the plot.    
                
        Notes    
        -----    
        This function generates a plot showing the total cross-section data and the best-fit curve.    
        If `plot_bg` is True, it will also plot the background function.    
        Can be used both after fitting (using fit_result) or before fitting (using model params).    
        """    
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig, ax = plt.subplots(figsize=(6, 4))    
            
        # Determine which results to use
        if stage is not None and hasattr(self, "fit_stages") and self.fit_stages:
            # Use specific stage results
            if stage < 1 or stage > len(self.fit_stages):
                raise ValueError(f"Stage {stage} not available. Available stages: 1-{len(self.fit_stages)}")
            
            # Get stage results
            stage_result = self.fit_stages[stage - 1]  # Convert to 0-indexed
            params = stage_result.params
            wavelength = np.linspace(1.0, 10.0, 1000)  # Adjust range as needed
            xs = self.cross_section(wavelength, **params)
            fit_label = f"Stage {stage} fit"
            
        elif hasattr(self, "fit_result") and self.fit_result is not None:    
            # Use final fit results    
            wavelength = self.fit_result.userkws["wl"]    
            params = self.fit_result.params    
            xs = self.cross_section(wavelength, **params)    
            fit_label = "Best fit"    
        else:    
            # Use model (no fit yet)    
            fit_label = "Model"    
            params = self.params    
            wavelength = np.linspace(1.0, 10.0, 1000)  # Adjust range as needed    
            xs = self.cross_section(wavelength, **params)    
            
        # Plot settings    
        color = kwargs.pop("color", "seagreen")    
        title = kwargs.pop("title", f"Total Cross-Section: {self.cross_section.name}")    
            
        # Plot cross-section    
        ax.plot(wavelength, xs, color=color, label=fit_label)    
        ax.set_ylabel("Cross-Section [barn]")    
        ax.set_xlabel("λ [Å]")    
        ax.set_title(title)    
            
        # Plot background if requested    
        if plot_bg and self.background:    
            bg = self.background.function(wl=wavelength, **params)    
            ax.plot(wavelength, bg, color="orange", linestyle="--", label="Background")    
            legend_labels = [fit_label, "Background"]    
        else:    
            legend_labels = [fit_label]    
            
        # Plot d-spacing lines if requested    
        if plot_dspace:    
            for phase in self.cross_section.phases_data:    
                try:    
                    hkls = self.cross_section.phases_data[phase].info.hklList()    
                except:    
                    continue    
                for hkl in hkls:    
                    hkl = hkl[:3]    
                    dspace = self.cross_section.phases_data[phase].info.dspacingFromHKL(*hkl)    
                    if dspace >= dspace_min:    
                        trans = ax.get_xaxis_transform()    
                        ax.axvline(dspace*2, lw=1, color="0.4", zorder=-1, ls=":")    
                        if len(self.cross_section.phases) > 1:    
                            ax.text(dspace*2, dspace_label_pos, f"{phase} {hkl}",     
                                    color="0.2", zorder=-1, fontsize=8, transform=trans,     
                                    rotation=90, va="top", ha="right")    
                        else:    
                            ax.text(dspace*2, dspace_label_pos, f"{hkl}",     
                                    color="0.2", zorder=-1, fontsize=8, transform=trans,     
                                    rotation=90, va="top", ha="right")    
            
        ax.legend(legend_labels, fontsize=9)    
        plt.tight_layout()    
        return ax

    def plot_stage_progression(self, param_name, ax=None, **kwargs):
        """
        Plot the progression of a parameter across fitting stages.

        Parameters
        ----------
        param_name : str
            The name of the parameter to plot (e.g., 'norm', 'thickness', 'b0').
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None, a new figure is created.
        **kwargs
            Additional keyword arguments for plotting (e.g., color, marker).

        Returns
        -------
        matplotlib.axes.Axes
            The axes containing the plot.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if not hasattr(self, 'fit_stages') or not self.fit_stages:
            raise ValueError("No stage results available. Run a multi-stage fit first.")

        if param_name not in self.params:
            raise ValueError(f"Parameter '{param_name}' not found. Available parameters: {list(self.params.keys())}")

        values = []
        stderrs = []
        stage_numbers = list(range(1, len(self.fit_stages) + 1))

        for stage_result in self.fit_stages:
            if param_name in stage_result.params:
                values.append(stage_result.params[param_name].value)
                stderrs.append(stage_result.params[param_name].stderr or 0)
            else:
                values.append(np.nan)
                stderrs.append(np.nan)

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))

        color = kwargs.pop("color", "seagreen")
        ax.errorbar(stage_numbers, values, yerr=stderrs, fmt="o-", color=color, **kwargs)
        ax.set_xlabel("Stage Number")
        ax.set_ylabel(f"{param_name}")
        ax.set_title(f"Progression of {param_name} Across Stages")
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        return ax

    def plot_chi2_progression(self, ax=None, **kwargs):
        """
        Plot the progression of reduced chi-squared across fitting stages.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None, a new figure is created.
        **kwargs
            Additional keyword arguments for plotting (e.g., color, marker).

        Returns
        -------
        matplotlib.axes.Axes
            The axes containing the plot.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if not hasattr(self, 'fit_stages') or not self.fit_stages:
            raise ValueError("No stage results available. Run a multi-stage fit first.")

        chi2_values = []
        stage_numbers = list(range(1, len(self.fit_stages) + 1))

        for stage_result in self.fit_stages:
            chi2_values.append(stage_result.redchi if hasattr(stage_result, 'redchi') else np.nan)

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))

        color = kwargs.pop("color", "seagreen")
        ax.plot(stage_numbers, chi2_values, "o-", color=color, **kwargs)
        ax.set_xlabel("Stage Number")
        ax.set_ylabel("Reduced χ²")
        ax.set_title("Reduced χ² Progression Across Stages")
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        return ax

    def get_stages_summary_table(self):
        """
        Get the stages summary table showing parameter progression through refinement stages.
        
        Returns
        -------
        pandas.DataFrame
            Multi-index DataFrame with parameters as rows and stages as columns.
            Each stage has columns for 'value', 'stderr', 'vary', and 'redchi'.
        """
        if not hasattr(self, "stages_summary"):
            raise ValueError("No stages summary available. Run fit with method='rietveld' first.")
        
        return self.stages_summary


    def interactive_plot(self, data=None, plot_bg=True, plot_dspace=False, 
                        dspace_min=1.0, dspace_label_pos=0.99, **kwargs):
        """
        Create an interactive plot with intuitive parameter controls using ipywidgets.

        Parameters
        ----------
        data : object, optional
            Data object to show alongside the model for comparison.
        plot_bg : bool, optional
            Whether to include the background in the plot, by default True.
        plot_dspace : bool, optional
            If True, plots 2*dspace lines and labels for materials with dspace >= dspace_min.
        dspace_min : float, optional
            Minimum dspace for plotting 2*dspace lines, by default 1.0.
        dspace_label_pos : float, optional
            Y-axis position for dspace labels, by default 0.99.
        kwargs : dict, optional
            Additional plot settings (e.g., color, marker size).

        Returns
        -------
        ipywidgets.VBox
            Container with interactive controls and plot.

        Notes
        -----
        Designed for models before fitting. Displays a warning if fit results exist.
        Provides real-time parameter exploration with sliders, float fields, and reset functionality.
        """
        # Check for fit results
        if hasattr(self, "fit_result") and self.fit_result is not None:
            print("Warning: interactive_plot is for models before fitting. Use plot() instead.")
            return

        # Store original parameters
        original_params = deepcopy(self.params)

        # Prepare data
        if data is not None:
            wavelength = data.table.wavelength
            data_values = data.table.trans
            err = data.table.err
        else:
            wavelength = np.linspace(1.0, 10.0, 1000)
            data_values = None
            err = None

        # Create output widget for plot
        plot_output = widgets.Output()

        # Dictionary for parameter widgets
        param_widgets = {}

        # Create parameter controls
        widget_list = []
        for param_name, param in self.params.items():
            # Parameter label
            label = widgets.Label(
                value=f"{param_name}:",
                layout={'width': '100px', 'padding': '5px'}
            )

            # Value slider
            if param.expr == "":
                slider = widgets.FloatSlider(
                    value=param.value,
                    min=param.min,
                    max=param.max,
                    # step=(param.max - param.min) / 2000,
                    readout=False,
                    disabled=not param.vary,
                    layout={'width': '200px'},
                    style={'description_width': '0px'}
                )
            else:
                slider = widgets.FloatSlider(
                    value=param.value,
                    min=0.001,  # For expressions, set a minimum to avoid zero division
                    max=1000,   # Arbitrary large max for expressions
                    step=(1000 - 0.001) / 200,
                    readout=False,
                    disabled=True,
                    layout={'width': '200px'},
                    style={'description_width': '0px'}
                )

            # Float text field
            float_text = widgets.FloatText(
                value=param.value,
                disabled=not param.vary,
                layout={'width': '80px'},
                style={'description_width': '0px'}
            )

            # Vary checkbox
            vary_widget = widgets.Checkbox(
                value=param.vary,
                description='Vary',
                layout={'width': '80px'},
                tooltip='Enable/disable parameter variation',
                style={'description_width': 'initial'}
            )

            # Store widgets
            param_widgets[param_name] = {'vary': vary_widget, 'float': float_text, 'slider': slider}

            # Create parameter row
            param_box = widgets.HBox([label, vary_widget, float_text, slider], layout={'padding': '2px'})
            widget_list.append(param_box)

            # Callbacks
            def make_update_callback(pname):
                def update_param(change):
                    # Sync slider and float text
                    if change['owner'] is param_widgets[pname]['slider']:
                        param_widgets[pname]['float'].value = change['new']
                    elif change['owner'] is param_widgets[pname]['float']:
                        param_widgets[pname]['slider'].value = change['new']
                    # Update model parameter
                    self.params[pname].value = param_widgets[pname]['slider'].value
                    self.params[pname].vary = param_widgets[pname]['vary'].value
                    # Enable/disable based on vary
                    if change['owner'] is param_widgets[pname]['vary']:
                        param_widgets[pname]['slider'].disabled = not change['new']
                        param_widgets[pname]['float'].disabled = not change['new']
                    # Update CrossSection with new parameters
                    param_kwargs = {pname: self.params[pname].value}
                    # Handle indexed parameters (e.g., ext_l1, a1) and non-indexed (e.g., α)
                    for param in self.params:
                        if param.endswith('1') or param in self.cross_section.materials:
                            param_kwargs[param] = self.params[param].value
                    self.cross_section(wavelength, **param_kwargs)
                    update_plot()
                return update_param

            slider.observe(make_update_callback(param_name), names='value')
            float_text.observe(make_update_callback(param_name), names='value')
            vary_widget.observe(make_update_callback(param_name), names='value')

        # Reset button
        reset_button = widgets.Button(
            description="Reset",
            button_style='info',
            tooltip='Reset parameters to original values',
            layout={'width': '100px'}
        )

        def reset_parameters(button):
            for param_name, original_param in original_params.items():
                self.params[param_name].value = original_param.value
                self.params[param_name].vary = original_param.vary
                param_widgets[param_name]['slider'].value = original_param.value
                param_widgets[param_name]['float'].value = original_param.value
                param_widgets[param_name]['vary'].value = original_param.vary
                param_widgets[param_name]['slider'].disabled = not original_param.vary
                param_widgets[param_name]['float'].disabled = not original_param.vary
            # Reset CrossSection with original parameters
            param_kwargs = {pname: original_params[pname].value for pname in original_params}
            self.cross_section(wavelength, **param_kwargs)
            update_plot()

        reset_button.on_click(reset_parameters)

        def update_plot():
            with plot_output:
                plot_output.clear_output(wait=True)
                model_values = self.eval(params=self.params, wl=wavelength)
                fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3.5, 1]}, figsize=(8, 6))

                # Plot settings
                color = kwargs.get("color", "teal")
                ecolor = kwargs.get("ecolor", "lightgray")
                title = kwargs.get("title", self.cross_section.name)
                ms = kwargs.get("ms", 2)

                # Plot data
                if data_values is not None:
                    residual = (data_values - model_values) / err
                    chi2 = np.sum(((data_values - model_values) / err) ** 2) / (len(data_values) - len(self.params))
                    ax0.errorbar(wavelength, data_values, err, marker="o", color=color, ms=ms, 
                                ecolor=ecolor, label="Data", zorder=1)
                    ax1.plot(wavelength, residual, color=color, linestyle='-', alpha=0.7)
                    chi2_text = f"χ²: {chi2:.2f}"
                else:
                    ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
                    chi2_text = "χ²: N/A"

                # Plot model
                ax0.plot(wavelength, model_values, color="navy", label="Model", linewidth=2, zorder=2)
                ax0.set_ylabel("Transmission", fontsize=10)
                ax0.set_title(title, fontsize=12, pad=10)

                ax1.set_ylabel("Residuals [1σ]", fontsize=10)
                ax1.set_xlabel("λ [Å]", fontsize=10)

                # Plot background
                if plot_bg and self.background:
                    self.background.plot(wl=wavelength, ax=ax0, params=self.params, **kwargs)
                    legend_labels = ["Model", "Background", "Data"] if data_values is not None else ["Model", "Background"]
                else:
                    legend_labels = ["Model", "Data"] if data_values is not None else ["Model"]

                # Legend
                ax0.legend(legend_labels, fontsize=9, loc='best', title=chi2_text, title_fontsize=9)

                # Plot d-spacing lines
                if plot_dspace:
                    for phase in self.cross_section.phases_data:
                        try:
                            hkls = self.cross_section.phases_data[phase].info.hklList()
                        except:
                            continue
                        for hkl in hkls:
                            hkl = hkl[:3]
                            dspace = self.cross_section.phases_data[phase].info.dspacingFromHKL(*hkl)
                            if dspace >= dspace_min:
                                ax0.axvline(dspace*2, lw=0.8, color="gray", ls=":", zorder=0)
                                trans = ax0.get_xaxis_transform()
                                label = f"{phase} {hkl}" if len(self.cross_section.phases) > 1 else f"{hkl}"
                                ax0.text(dspace*2, dspace_label_pos, label, color="darkgray", fontsize=8, 
                                        transform=trans, rotation=90, va="top", ha="right")

                plt.subplots_adjust(hspace=0.05)
                plt.tight_layout()
                plt.show()

        # Layout
        controls_box = widgets.VBox(
            [widgets.HTML("<h4 style='margin: 5px;'>Parameter Controls</h4>"), reset_button] + widget_list,
            layout={'padding': '10px', 'border': '1px solid lightgray', 'width': '350px'}
        )
        main_box = widgets.HBox([controls_box, plot_output])

        # Initial plot
        update_plot()
        return main_box

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

        return self