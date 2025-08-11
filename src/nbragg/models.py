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
from typing import List, Optional
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
            It True, allows the lattice parameters of the material to be varied 
        vary_extinction: bool, optional
            It True, allows the extinction parameters of the material to be varied (requires the CrysExtn plugin to be installed)
        kwargs : dict, optional
            Additional keyword arguments for model and background parameters.

        Notes
        -----
        This model calculates the transmission function as a combination of 
        cross-section, response function, and background.
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

        if params!=None:
            self.params = params.copy()
        else:
            self.params = lmfit.Parameters()
        if "thickness" not in self.params and "norm" not in self.params:
            self.params += self._make_basic_params()
        if "temp" not in self.params:
            self.params += self._make_temperature_params() # add temperature params to fit
        if vary_weights is not None:
            self.params += self._make_weight_params(vary=vary_weights)
        if vary_tof is not None:
            self.params += self._make_tof_params(vary=vary_tof,**kwargs)
        if vary_lattice is not None:
            self.params += self._make_lattice_params(vary=vary_lattice)
        if vary_extinction is not None:
            self.params += self._make_extinction_params(vary=vary_extinction)


        self.response = None
        if vary_response is not None:
            self.response = Response(kind=response,vary=vary_response)
            if list(self.response.params.keys())[0] in self.params:
                for param_name in self.params.keys():
                    self.params[param_name].vary = vary_response 
            else:
                self.params += self.response.params


        self.background = None
        if vary_background is not None:
            self.background = Background(kind=background,vary=vary_background)
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
            method: str = "leastsq",
            xtol: float = None, ftol: float = None, gtol: float = None,
            verbose: bool = False,
            progress_bar: bool = True,
            param_groups: Optional[List[List[str]]] = None,
            **kwargs):
        """
        Fit the model to data.

        This method supports both:
        - **Standard single-stage fitting** (default)
        - **Rietveld-style staged refinement** (`method="rietveld"`)

        Parameters
        ----------
        data : pandas.DataFrame or Data or array-like
            The input data.  
            - For `pandas.DataFrame` or `Data`: must have columns `"wavelength"`, `"trans"`, and `"err"`.
            - For array-like: will be passed directly to `lmfit.Model.fit`.
        params : lmfit.Parameters, optional
            Parameters to use for fitting. If None, uses the model's default parameters.
        wlmin, wlmax : float, optional
            Minimum and maximum wavelength for fitting (ignored for array-like input).
        method : str, optional
            Fitting method.  
            - `"leastsq"` (default) or any method supported by `lmfit`.
            - `"rietveld"` will run staged refinement via `_rietveld_fit`.
        xtol, ftol, gtol : float, optional
            Convergence tolerances (passed to `lmfit`).
        verbose : bool, optional
            If True, prints detailed fitting information.
        progress_bar : bool, optional
            If True, shows a progress bar for fitting:
            - For `"rietveld"`: shows stage name and reduced chi² per stage.
            - For regular fits: shows overall fit progress.
        param_groups : list, dict, or None, optional
            Used only for `"rietveld"`. Groups of parameters to fit in each stage.
            See `_rietveld_fit` docstring for details.
        **kwargs
            Additional keyword arguments passed to `lmfit.Model.fit`.

        Returns
        -------
        lmfit.model.ModelResult
            The fit result object, with extra methods:
            - `.plot()` — plot the fit result.
            - `.plot_total_xs()`, `.plot_stage_progression()`, `.plot_chi2_progression()` for advanced diagnostics.
            - `.stages_summary` (for `"rietveld"`).

        Examples
        --------
        **Basic fit:**
        ```python
        result = model.fit(data_df, wlmin=1.0, wlmax=5.0)
        result.plot()
        ```

        **Rietveld-style staged refinement:**
        ```python
        param_groups = {
            "Norm/Thick": ["norm", "thickness"],
            "Background": ["b0", "b1"],
            "Extinction": ["ext_l", "ext_Gg"]
        }
        result, summary = model.fit(
            data_df, method="rietveld",
            param_groups=param_groups,
            progress_bar=True,
            return_all_results=True
        )
        print(summary)
        ```

        Notes
        -----
        - `"rietveld"` mode is a **staged refinement** approach where parameters are refined in
        groups, improving stability for complex models.
        - Progress bars use `tqdm` and will automatically adapt to Jupyter or terminal output.
        """
        # Route to Rietveld if requested
        if method == "rietveld":
            return self._rietveld_fit(
                data, params, wlmin, wlmax,
                verbose=verbose,
                progress_bar=progress_bar,
                param_groups=param_groups,
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


    def _rietveld_fit(self, data, params: "lmfit.Parameters" = None, wlmin:float=1, wlmax:float=8,
                    verbose=False, progress_bar=True,
                    param_groups=None,
                    **kwargs):
        """ Perform Rietveld-style staged fitting.
        This method allows for staged fitting of parameters, where each stage
        refines a specific group of parameters while keeping others fixed.
        Parameters
        ----------
        data : pandas.DataFrame or Data
            The input data containing wavelength and transmission values.
        params : lmfit.Parameters, optional
            Initial parameters for the fit. If None, uses the model's default parameters.       
        wlmin : float, optional default=1
            Minimum wavelength for fitting.
        wlmax : float, optional default=8
            Maximum wavelength for fitting.
        verbose : bool, optional
            If True, prints detailed information about each fitting stage.
        progress_bar : bool, optional
            If True, shows a progress bar for each fitting stage.
        param_groups : list, dict, or None, optional - only used for Rietveld fitting
            Groups of parameters to fit in each stage. Can be:
            - List of lists: [["norm", "thickness"], ["background", "extinction"]]
            - List of strings/lists: ["basic", ["b0", "ext_l2"]]
            - Dict with stage names: {"stage1": ["norm"], "stage2": ["background"]}
            If None, uses predefined groups.
        kwargs : dict, optional
            Additional keyword arguments for the fit method, such as weights, method, etc.

        Returns
        -------
        fit_result : lmfit.ModelResult
            The final fit result after all stages.

        fit_result.stages_summary : pandas.DataFrame 
            Summary of each fitting stage, including parameter values and reduced chi-squared.

        Notes
        -----
        This method is designed for Rietveld-style fitting, where parameters are fitted in stages.
        It allows for flexible grouping of parameters, enabling users to refine specific aspects of the model iteratively.
        """
        from copy import deepcopy
        import sys
        import warnings
        try:
            from tqdm.notebook import tqdm
        except ImportError:
            from tqdm.auto import tqdm
        import pickle

        # User-friendly group name mapping
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p) ],
            "lattice": [p for p in self.params if p in ["a", "b", "c"]],
            "extinction": [p for p in self.params if p.startswith("ext_")],
            "orientation": [p for p in self.params if p.startswith("θ") or p.startswith("ϕ") or p.startswith("η")],
            "mosaicity": [p for p in self.params if p.startswith("η")],
            "temperature": [p for p in ["temp"] if p in self.params],
        }

        def resolve_single_param_or_group(item):
            """Resolve a single parameter name or group name to a list of parameters."""
            if item in group_map:
                # It's a predefined group
                resolved = group_map[item]
                if verbose:
                    print(f"  Resolved group '{item}' to: {resolved}")
                return resolved
            elif item in self.params:
                # It's a single parameter name
                if verbose:
                    print(f"  Found parameter: {item}")
                return [item]
            else:
                # Check if it matches any parameters with wildcards
                matching_params = [p for p in self.params.keys() if fnmatch.fnmatch(p, item)]
                if matching_params:
                    if verbose:
                        print(f"  Pattern '{item}' matched: {matching_params}")
                    return matching_params
                else:
                    warnings.warn(f"Unknown parameter or group: '{item}'. Available parameters: {list(self.params.keys())}")
                    return []

        def resolve_group(entry):
            """Resolve a group entry (string, list, or nested structure) to a flat list of parameters."""
            if isinstance(entry, str):
                return resolve_single_param_or_group(entry)
            elif isinstance(entry, list):
                # Flatten nested lists
                resolved = []
                for item in entry:
                    if isinstance(item, str):
                        resolved.extend(resolve_single_param_or_group(item))
                    elif isinstance(item, list):
                        # Handle nested lists recursively
                        resolved.extend(resolve_group(item))
                    else:
                        warnings.warn(f"Unexpected item type in group: {type(item)} - {item}")
                return resolved
            else:
                warnings.warn(f"Invalid param_groups entry type: {type(entry)} - {entry}")
                return []

        # Handle different input formats for param_groups
        stage_names = []
        if param_groups is None:
            # Default groups
            param_groups = [
                "basic", "background", "tof", "response",
                "weights", "lattice", "extinction", "orientation","mosaicity", "temperature",
            ]
            resolved_param_groups = [resolve_group(g) for g in param_groups]
            stage_names = [f"Stage_{i+1}" for i in range(len(param_groups))]
            
        elif isinstance(param_groups, dict):
            # Dictionary format: {"stage_name": ["param1", "param2"], ...}
            stage_names = list(param_groups.keys())
            resolved_param_groups = [resolve_group(param_groups[stage]) for stage in stage_names]
            if verbose:
                print(f"Using custom stage names: {stage_names}")
                
        elif isinstance(param_groups, list):
            # List format: [["param1", "param2"], ["param3"], ...]
            resolved_param_groups = [resolve_group(g) for g in param_groups]
            stage_names = [f"Stage_{i+1}" for i in range(len(param_groups))]
            
        else:
            raise ValueError("param_groups must be None, a list, or a dictionary")

        # Remove empty groups
        valid_groups = []
        valid_names = []
        for i, group in enumerate(resolved_param_groups):
            if group:  # Only keep non-empty groups
                valid_groups.append(group)
                valid_names.append(stage_names[i])
            elif verbose:
                print(f"Skipping empty group: {stage_names[i]}")
        
        resolved_param_groups = valid_groups
        stage_names = valid_names

        if not resolved_param_groups:
            raise ValueError("No valid parameter groups found. Check your parameter names.")

        if verbose:
            print(f"\nFitting stages:")
            for i, (name, group) in enumerate(zip(stage_names, resolved_param_groups)):
                print(f"  {name}: {group}")

        # Store the resolved param groups for the summary table
        self._stage_param_groups = resolved_param_groups
        self._stage_names = stage_names

        # Prepare data
        if isinstance(data, pandas.DataFrame):
            data = data.query(f"{wlmin} < wavelength < {wlmax}")
            wavelengths = data["wavelength"].values
            trans = data["trans"].values
            weights = kwargs.get("weights", 1. / data["err"].values)
        elif isinstance(data, Data):
            data = data.table.query(f"{wlmin} < wavelength < {wlmax}")
            wavelengths = data["wavelength"].values
            trans = data["trans"].values
            weights = kwargs.get("weights", 1. / data["err"].values)
        else:
            raise ValueError("Rietveld fitting requires wavelength-based input data.")

        params = deepcopy(params or self.params)

        # Use tqdm.notebook for Jupyter environments
        try:
            from tqdm.notebook import tqdm as notebook_tqdm
            # Check if we're in a Jupyter environment
            if 'ipykernel' in sys.modules:
                iterator = notebook_tqdm(
                    zip(stage_names, resolved_param_groups), 
                    desc="Rietveld Fit", 
                    disable=not progress_bar,
                    total=len(stage_names)
                )
            else:
                iterator = tqdm(
                    zip(stage_names, resolved_param_groups), 
                    desc="Rietveld Fit", 
                    disable=not progress_bar,
                    total=len(stage_names)
                )
        except ImportError:
            iterator = tqdm(
                zip(stage_names, resolved_param_groups), 
                desc="Rietveld Fit", 
                disable=not progress_bar,
                total=len(stage_names)
            )
        
        stage_results = []  # Store results
        stage_summaries = []  # For DataFrame summary

        def extract_pickleable_attributes(fit_result):
            """Extract only pickleable attributes from fit_result"""
            # List of commonly pickleable attributes from lmfit ModelResult
            safe_attrs = [
                'params', 'success', 'residual', 'chisqr', 'redchi', 'aic', 'bic',
                'nvarys', 'ndata', 'nfev', 'message', 'lmdif_message', 'cov_x',
                'method', 'flatchain', 'errorbars', 'ci_out'
            ]
            
            # Create a simple object to hold the results
            class PickleableResult:
                def __init__(self):
                    pass
            
            result = PickleableResult()
            
            for attr in safe_attrs:
                if hasattr(fit_result, attr):
                    try:
                        value = getattr(fit_result, attr)
                        # Test if it's pickleable
                        pickle.dumps(value)
                        setattr(result, attr, value)
                    except (TypeError, ValueError, AttributeError):
                        # Skip non-pickleable attributes
                        if verbose:
                            print(f"Skipping non-pickleable attribute: {attr}")
                        continue
            
            return result

        for stage_name, group in iterator:
            stage_num = len(stage_results) + 1

            if verbose:
                print(f"\n{stage_name}: Fitting {group}")

            # Freeze all parameters
            for p in params.values():
                p.vary = False

            # Unfreeze current group
            unfrozen_count = 0
            for name in group:
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

            # Fit this stage
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

            # Extract only pickleable parts
            stripped_result = extract_pickleable_attributes(fit_result)

            # Store results
            stage_results.append(stripped_result)

            # Build summary row
            summary = {
                "stage": stage_num,
                "stage_name": stage_name,
                "fitted_params": group,
                "redchi": fit_result.redchi
            }
            for name, par in fit_result.params.items():
                summary[f"{name}_value"] = par.value
                summary[f"{name}_stderr"] = par.stderr
                summary[f"{name}_vary"] = par.vary
            stage_summaries.append(summary)

            # Update tqdm display
            iterator.set_description(f"Stage {stage_num}/{len(stage_names)}")
            iterator.set_postfix({"stage": stage_name, "reduced χ²": f"{fit_result.redchi:.4g}"})

            # Update for next stage
            params = fit_result.params

            if verbose:
                print(f"  {stage_name} completed. χ²/dof = {fit_result.redchi:.4f}")


        if not stage_results:
            raise RuntimeError("No successful fitting stages completed")

        # Final
        self.fit_result = fit_result
        self.fit_stages = stage_results
        
        # Call the updated summary table method with stage names
        self.stages_summary = self._create_stages_summary_table_enhanced(stage_results, resolved_param_groups, stage_names)

        # Attach plotting
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

    def _create_stages_summary_table_enhanced(self, stage_results, resolved_param_groups, stage_names=None, color=True):
        import pandas as pd
        import numpy as np

        # --- Build the DataFrame ---
        all_param_names = list(stage_results[-1].params.keys())
        stage_data = {}
        if stage_names is None:
            stage_names = [f"Stage_{i+1}" for i in range(len(stage_results))]

        for stage_idx, stage_result in enumerate(stage_results):
            stage_col = stage_names[stage_idx] if stage_idx < len(stage_names) else f"Stage_{stage_idx + 1}"
            stage_data[stage_col] = {'value': {}, 'stderr': {}, 'vary': {}}
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

        # 1) Highlight vary=True cells (light blue)
        vary_cols = [col for col in df.columns if col[1] == 'vary']
        def highlight_vary(s):
            return ['background-color: lightblue' if v is True else '' for v in s]
        for col in vary_cols:
            styler = styler.apply(highlight_vary, subset=[col], axis=0)

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
        if show_groups:
            print("Available parameter groups:")
            print("=" * 30)

            group_map = {
                "basic": ["norm", "thickness"],
                "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
                "tof": [p for p in ["L0", "t0"] if p in self.params],
                "response": [p for p in self.params if self.response and p in self.response.params],
                "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
                "lattice": [p for p in self.params if p in ["a", "b", "c"]],
                "extinction": [p for p in self.params if p.startswith("ext_")],
                "orientation": [p for p in ["phi", "theta", "eta"] if p in self.params],
                "temperature": [p for p in ["temp"] if p in self.params],
            }
            
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
        This function generates a plot showing the total cross-section data,     
        the best-fit curve, and residuals. If `plot_bg` is True, it will also     
        plot the background function.    
        """    
        # Determine which parameters to use
        if stage is not None and hasattr(self, "fit_stages") and self.fit_stages:
            # Use specific stage results
            if stage < 1 or stage > len(self.fit_stages):
                raise ValueError(f"Stage {stage} not available. Available stages: 1-{len(self.fit_stages)}")
            
            stage_result = self.fit_stages[stage - 1]  # Convert to 0-indexed
            params = stage_result.params
            chi2 = stage_result.redchi if hasattr(stage_result, 'redchi') else np.nan
            title_suffix = f" (Stage {stage})"
            
            # We still need the original wavelength data
            if not hasattr(self, "fit_result") or self.fit_result is None:
                raise ValueError("Cannot plot stage cross-section without original fit data")
            wavelength = self.fit_result.userkws["wl"]
            data_values = self.fit_result.data
            
        elif hasattr(self, "fit_result") and self.fit_result is not None:
            # Use final fit results
            params = self.fit_result.params
            chi2 = self.fit_result.redchi
            wavelength = self.fit_result.userkws["wl"]
            data_values = self.fit_result.data
            title_suffix = ""
        else:
            raise ValueError("Cannot plot cross-section without fit results")
        
        # Prepare data for cross-section calculation    
        if "k" in params:    
            k = params['k'].value    
        else:    
            k = 1.    
        if plot_bg and self.background:    
            bg = self.background.function(wavelength, **params)    
        else:    
            bg = 0.    
        norm = params['norm'].value    
        n = self.atomic_density    
        thickness = params['thickness'].value    

        # Calculate cross-section data    
        data_xs = -1. / n / thickness * np.log((data_values - k * bg) / norm / (1. - bg))    

        fig, ax = plt.subplots(2, 1, sharex=True, height_ratios=[3.5, 1], figsize=(6, 5))    

            
        # Calculate best fit and residuals for cross-section    
        xs = self.cross_section(wavelength, **params)  # You'll need to implement this method    

        if self.response != None:    
            response = self.response.function(**params)    
            best_fit = convolve1d(xs, response, 0)    
        else:
            best_fit = xs
        residual = data_xs - best_fit    

        # Plot styling    
        color = kwargs.pop("color", "crimson")    
        title = kwargs.pop("title", self.cross_section.name)    
        ecolor = kwargs.pop("ecolor", "0.8")    
        ms = kwargs.pop("ms", 2)    

        # Top subplot: Data and fit    
        ax[0].errorbar(wavelength, data_xs,     
                    # yerr=1./self.fit_result.weights,  # Assuming similar error handling    
                    marker="o",     
                    color=color,     
                    ms=ms,     
                    zorder=-1,     
                    ecolor=ecolor,     
                    label="Cross-section data")    
            
        ax[0].plot(wavelength, best_fit, color="0.4", label="Best fit")    
        ax[0].plot(wavelength, xs, color="0.2", label="total xs")    
        ax[0].set_ylabel("Total Cross Section [barn/sr]")    
        ax[0].set_title(title)    

        # Bottom subplot: Residuals    
        ax[1].plot(wavelength, residual, color=color)    
        ax[1].set_ylabel("Residuals [1σ]")    
        ax[1].set_xlabel("λ [Å]")    

        # Background plotting (if enabled)    
        if plot_bg and self.background:    
            self.background.plot(wl=wavelength, ax=ax[0], params=params, **kwargs)    
            ax[0].legend(["Cross-section data", "Background", "Total cross-section","Best fit"][::-1],     
                        fontsize=9,     
                        reverse=True,     
                        title=f"χ$^2$: {chi2:.2f}" if not np.isnan(chi2) else "χ$^2$: N/A")    
        else:    
            ax[0].legend(["Cross-section data","Total cross-section", "Best fit"][::-1],     
                        fontsize=9,     
                        reverse=True,     
                        title=f"χ$^2$: {chi2:.2f}" if not np.isnan(chi2) else "χ$^2$: N/A")    

        # d-spacing plot (if enabled)    
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
                            
                        # Label d-spacing lines    
                        if len(self.cross_section.phases) > 1:    
                            ax[0].text(dspace*2, dspace_label_pos,     
                                    f"{phase} {hkl}",     
                                    color="0.2",     
                                    zorder=-1,     
                                    fontsize=8,     
                                    transform=trans,     
                                    rotation=90,     
                                    va="top",     
                                    ha="right")    
                        else:    
                            ax[0].text(dspace*2, dspace_label_pos,     
                                    f"{hkl}",     
                                    color="0.2",     
                                    zorder=-1,     
                                    fontsize=8,     
                                    transform=trans,     
                                    rotation=90,     
                                    va="top",     
                                    ha="right")    

        plt.subplots_adjust(hspace=0.05)    
        return ax

    def plot_stage_progression(self, stages: list = None, **kwargs):
        """
        Plot the progression of Rietveld refinement stages showing how the fit improves.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if not hasattr(self, "fit_stages") or not self.fit_stages:
            raise ValueError("No Rietveld stages available. Run fit with method='rietveld' first.")

        if stages is None:
            stages = list(range(1, len(self.fit_stages) + 1))

        # Original data
        if hasattr(self, "fit_result") and self.fit_result is not None:
            wavelength = self.fit_result.userkws["wl"]
            data_values = self.fit_result.data
            err = 1. / self.fit_result.weights
        else:
            raise ValueError("Cannot plot stage progression without original fit data")

        fig, ax = plt.subplots(figsize=(6, 4))

        # Match style: light gray points for data
        ax.errorbar(wavelength, data_values, err,
                    marker="o", color="0.6", ms=2, alpha=0.7, zorder=-1,
                    ecolor="0.85", label="Data")

        # Use consistent style palette
        colors = plt.cm.plasma(np.linspace(0, 0.85, len(stages)))

        for i, stage in enumerate(stages):
            if stage < 1 or stage > len(self.fit_stages):
                continue

            stage_result = self.fit_stages[stage - 1]
            params = stage_result.params
            best_fit = self.eval(params=params, wl=wavelength)
            chi2 = getattr(stage_result, "redchi", np.nan)

            # Get stage name if available
            stage_name = f"Stage {stage}"
            if hasattr(self, "stages_summary"):
                stage_col = f"Stage_{stage}"
                if (stage_col, "vary") in self.stages_summary.columns:
                    varied_params = self.stages_summary.loc[
                        self.stages_summary[(stage_col, "vary")] == True
                    ].index.tolist()
                    varied_params = [p for p in varied_params if p != "redchi"]
                    if varied_params:
                        stage_name = ", ".join(varied_params[:2]) + (
                            f" +{len(varied_params)-2}" if len(varied_params) > 2 else ""
                        )

            ax.plot(wavelength, best_fit,
                    color=colors[i], lw=1.2 + 0.4 * i,
                    alpha=0.8,
                    label=f"{stage_name} (χ²={chi2:.3f})" if not np.isnan(chi2) else stage_name)

        ax.set_xlabel("λ [Å]")
        ax.set_ylabel("Transmission")
        ax.set_title("Rietveld Refinement Stage Progression")
        ax.legend(fontsize=8, frameon=False)

        plt.tight_layout()
        return ax


    def plot_chi2_progression(self, **kwargs):
        """
        Plot the χ² progression through Rietveld stages with stage names on x-axis.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if not hasattr(self, "fit_stages") or not self.fit_stages:
            raise ValueError("No Rietveld stages available. Run fit with method='rietveld' first.")

        stages = list(range(1, len(self.fit_stages) + 1))
        chi2_values = []
        stage_labels = []

        for stage in stages:
            stage_result = self.fit_stages[stage - 1]
            chi2 = getattr(stage_result, "redchi", np.nan)
            chi2_values.append(chi2)

            label = f"Stage {stage}"
            if hasattr(self, "stages_summary"):
                stage_col = f"Stage_{stage}"
                if (stage_col, "vary") in self.stages_summary.columns:
                    varied_params = self.stages_summary.loc[
                        self.stages_summary[(stage_col, "vary")] == True
                    ].index.tolist()
                    varied_params = [p for p in varied_params if p != "redchi"]
                    if varied_params:
                        label = ", ".join(varied_params[:2]) + (
                            f" +{len(varied_params)-2}" if len(varied_params) > 2 else ""
                        )
            stage_labels.append(label)

        fig, ax = plt.subplots(figsize=(6, 3.5))

        ax.plot(stages, chi2_values, marker="o", lw=2, color="seagreen")

        # Annotate each point
        for stage, chi2 in zip(stages, chi2_values):
            if not np.isnan(chi2):
                ax.annotate(f"{chi2:.3f}", (stage, chi2),
                            textcoords="offset points", xytext=(0, 8),
                            ha="center", fontsize=8)

        ax.set_xlabel("Refinement Stage")
        ax.set_ylabel("Reduced χ²")
        ax.set_title("Rietveld χ² Progression")

        # Stage names at bottom
        ax.set_xticks(stages)
        ax.set_xticklabels(stage_labels, rotation=30, ha="right", fontsize=8)

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
        
    def _make_orientation_params(self, vary: bool = False):
        """
        Create orientation for the model.

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
        for i, material in enumerate(self._materials):
            mos = self._materials[material].get("mos", None)
            if mos: 
                param_name = f"η{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=mos, min=0.001, max=50, vary=vary)
                
                theta = self._materials[material].get("theta", 0.)
                param_name = f"θ{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=theta, min=0., max=180, vary=vary)
                
                phi = self._materials[material].get("phi", 0.)
                param_name = f"ϕ{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=phi, min=0, max=360, vary=vary)
        return params

    def _make_temperature_params(self, vary: bool = False):
        """
        Create temperature for the model.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The temperature-related parameters.
        """
        params = lmfit.Parameters()
        for i, material in enumerate(self._materials):
            temp = self._materials[material].get("temp", None)
            if temp: 
                param_name = "temp"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=temp, min=77., max=1000, vary=vary)
        return params

    def _make_basic_params(self, vary: bool = True):
        """
        Create basic params for the model.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The basic parameters of thickness and normalization.
        """
        params = lmfit.Parameters()
        param_name = "thickness"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=1., min=0., max=5., vary=vary)
        
        param_name = "norm"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=1., min=0.1, max=10., vary=vary)
        return params

    def _make_tof_params(self, vary: bool = False, t0: float = 0., L0: float = 1.):
        """
        Create time-of-flight (TOF) parameters for the model.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow these parameters to vary during fitting, by default False.
        t0 : float, optional
            Initial time offset parameter, by default 0.
        L0 : float, optional
            Initial flight path distance scale parameter, by default 1.

        Returns
        -------
        lmfit.Parameters
            The TOF-related parameters.
        """
        params = lmfit.Parameters()
        param_name = "L0"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=L0, min=0.5, max=1.5, vary=vary)
        
        param_name = "t0"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=t0, min=-5e-6, max=5e6, vary=vary)
        return params

    def _make_weight_params(self, vary: bool = False):
        """
        Create lmfit parameters based on initial isotope weights.

        Parameters
        ----------
        vary : bool, optional
            Whether to allow weights to vary during fitting, by default False.

        Returns
        -------
        lmfit.Parameters
            The normalized weight parameters for the model.
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
    
    def _make_lattice_params(self, vary: bool = False):
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

    
    def _make_extinction_params(self, vary: bool = False):
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