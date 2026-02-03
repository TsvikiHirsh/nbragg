"""
Fitting-related methods for TransmissionModel.

This module contains the FittingMixin class which provides fitting functionality
for neutron transmission models, including single-stage, multi-stage (Rietveld),
and grouped/parallel fitting capabilities.
"""

import lmfit
import numpy as np
import pandas
import warnings
import re
import fnmatch
import pickle
from copy import deepcopy
from typing import Dict, List, Optional, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from nbragg.data import Data
    from nbragg.cross_section import CrossSection
    from nbragg.response import Response, Background

from nbragg.data import Data
from nbragg.grouped_fit import (
    GroupedFitResult,
    _fit_single_group_worker,
    _reconstruct_result_from_dict,
    _add_save_method_to_result
)


class FittingMixin:
    """
    Mixin class providing fitting functionality for TransmissionModel.

    This mixin provides:
    - Single-stage least-squares fitting
    - Multi-stage Rietveld and staged refinement
    - Parallel fitting for grouped/spatially-resolved data
    - Support for multiple parallelization backends (loky, threading, sequential)
    """

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
        n_jobs : int, optional
            Number of parallel jobs for grouped data fitting (default: 10).
            Only applies when fitting grouped data. Set to 1 for sequential fitting.
            Set to -1 to use all available CPU cores.
        backend : str, optional
            Parallelization backend for grouped data fitting (default: "loky").
            Options:
            - "loky": True multiprocessing with model reconstruction in each worker.
              Provides full CPU parallelism but has overhead from recreating NCrystal
              objects in each process.
            - "threading": Threading-based parallelism (limited by Python's GIL).
              Lower overhead but limited speedup for CPU-bound tasks.
            - "sequential": No parallelization. Useful for debugging.
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

        # For grouped data with parallel fitting
        >>> grouped_result = model.fit(grouped_data, n_jobs=10)
        >>> result_0_0 = grouped_result[(0, 0)]
        >>> grouped_result.plot_parameter_map("thickness")
        """
        # Check if data is grouped and route to parallel fitting
        if hasattr(data, 'is_grouped') and data.is_grouped:
            n_jobs = kwargs.pop('n_jobs', 10)
            backend = kwargs.pop('backend', 'loky')
            return self._fit_grouped(
                data, params, wlmin, wlmax,
                method=method,
                xtol=xtol, ftol=ftol, gtol=gtol,
                verbose=verbose,
                progress_bar=progress_bar,
                stages=stages,
                n_jobs=n_jobs,
                backend=backend,
                **kwargs
            )

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

        # Check if lattice parameters are active - if so, increase finite-difference step.
        # NCrystal requires lattice parameter changes >= ~1e-4 A to produce
        # different cross-sections, but lmfit's default Jacobian step is ~1.5e-8
        lattice_pattern = re.compile(r'^[abc]\d*$')
        fit_params = params or self.params
        active_lattice = any(
            lattice_pattern.match(name) and fit_params[name].vary
            for name in fit_params
        )
        if active_lattice:
            if method in ("leastsq",):
                fit_kws.setdefault("epsfcn", 1e-4)
            elif method in ("least_squares", "least-squares"):
                fit_kws.setdefault("diff_step", 1e-4)

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

        # Add save() method to the result
        return _add_save_method_to_result(fit_result)

    def _get_stage_parameters(self, stage_def: Union[str, List[str]]) -> List[str]:
        """Helper method to get parameters associated with a stage definition."""
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
            "lattice": [p for p in self.params if re.compile(r'^[abc]\d*$').match(p)],
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
        import sys

        try:
            from tqdm.notebook import tqdm
        except ImportError:
            from tqdm.auto import tqdm

        if method not in ["rietveld", "staged"]:
            raise ValueError(f"Invalid multi-stage method: {method}. Use 'rietveld' or 'staged'.")

        # User-friendly group name mapping
        group_map = {
            "basic": ["norm", "thickness"],
            "background": [p for p in self.params if re.compile(r"(b|bg)\d+").match(p) or p.startswith("b_")],
            "tof": [p for p in ["L0", "t0"] if p in self.params],
            "response": [p for p in self.params if self.response and p in self.response.params],
            "weights": [p for p in self.params if re.compile(r"p\d+").match(p)],
            "lattice": [p for p in self.params if re.compile(r'^[abc]\d*$').match(p)],
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
                # But respect the user's vary setting from self.params
                unfrozen_count = 0
                for name in cumulative_params:
                    if name in params:
                        # Only set vary=True if the parameter's original setting allows it
                        if name in self.params and self.params[name].vary:
                            params[name].vary = True
                            unfrozen_count += 1
                            if verbose and (name in group or not group):
                                print(f"  New parameter: {name}")
                            elif verbose:
                                print(f"  Continuing: {name}")
                        elif verbose:
                            print(f"  Skipping {name} (vary=False set by user)")
                    else:
                        if name in group or not group:  # Only warn for new parameters
                            warnings.warn(f"Parameter '{name}' not found in params")

                if verbose:
                    print(f"  Total active parameters: {unfrozen_count}")

            elif method == "staged":
                # Staged: only current group parameters vary
                # But respect the user's vary setting from self.params
                for p in params.values():
                    p.vary = False

                unfrozen_count = 0
                active_params = group if group else [p for p in self.params if self.params[p].vary]
                for name in active_params:
                    if name in params:
                        # Only set vary=True if the parameter's original setting allows it
                        if name in self.params and self.params[name].vary:
                            params[name].vary = True
                            unfrozen_count += 1
                            if verbose:
                                print(f"  Unfrozen: {name}")
                        elif verbose:
                            print(f"  Skipping {name} (vary=False set by user)")
                    else:
                        warnings.warn(f"Parameter '{name}' not found in params")

            if unfrozen_count == 0:
                warnings.warn(f"No parameters were unfrozen in {stage_name}. Skipping this stage.")
                continue

            # Check if lattice parameters are active - if so, increase epsfcn
            # NCrystal requires lattice parameter changes >= ~1e-4 A to produce
            # different cross-sections, but lmfit's default Jacobian step is ~1.5e-8
            lattice_pattern = re.compile(r'^[abc]\d*$')
            active_lattice = any(
                lattice_pattern.match(name) and params[name].vary
                for name in params if name in params
            )
            stage_kwargs = dict(kwargs)
            if active_lattice:
                fit_kws = stage_kwargs.pop("fit_kws", {})
                fit_kws.setdefault("epsfcn", 1e-4)
                stage_kwargs["fit_kws"] = fit_kws

            # Perform fitting
            try:
                fit_result = super().fit(
                    trans,
                    params=params,
                    wl=wavelengths,
                    weights=weights,
                    method="leastsq",
                    **stage_kwargs
                )
            except Exception as e:
                if verbose:
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
        fit_result.fit_stages = stage_results  # Also store on result for grouped data access
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

        # Add save() method to the result
        return _add_save_method_to_result(fit_result)

    def _fit_grouped(self, data, params=None, wlmin: float = 1., wlmax: float = 6.,
                     method: str = "rietveld",
                     xtol: float = None, ftol: float = None, gtol: float = None,
                     verbose: bool = False,
                     progress_bar: bool = True,
                     stages: Optional[Union[str, Dict[str, Union[str, List[str]]]]] = None,
                     n_jobs: int = 10,
                     backend: str = "loky",
                     **kwargs):
        """
        Fit model to grouped data in parallel.

        Parameters:
        -----------
        data : Data
            Grouped data object with is_grouped=True.
        params : lmfit.Parameters, optional
            Parameters to use for fitting.
        wlmin, wlmax : float
            Wavelength range for fitting.
        method : str
            Fitting method: "least-squares", "rietveld", or "staged".
        xtol, ftol, gtol : float, optional
            Convergence tolerances.
        verbose : bool
            Show progress for individual fits.
        progress_bar : bool
            Show overall progress bar.
        stages : str or dict, optional
            Fitting stages configuration.
        n_jobs : int
            Number of parallel jobs (default: 10).
        backend : str
            Parallelization backend: "loky" (true multiprocessing, default),
            "threading" (GIL-limited but works with shared objects),
            or "sequential" (no parallelization, for debugging).
        **kwargs
            Additional arguments passed to fit.

        Returns:
        --------
        GroupedFitResult
            Container with fit results for each group.
        """
        if backend == "loky":
            return self._fit_grouped_loky(data, params, wlmin, wlmax, method,
                                          xtol, ftol, gtol, verbose, progress_bar,
                                          stages, n_jobs, **kwargs)
        elif backend == "threading":
            if n_jobs > 4:
                print(f"      Consider n_jobs=4 or less for better performance with threading.")
            return self._fit_grouped_threading(data, params, wlmin, wlmax, method,
                                               xtol, ftol, gtol, verbose, progress_bar,
                                               stages, n_jobs, **kwargs)
        elif backend == "sequential":
            return self._fit_grouped_sequential(data, params, wlmin, wlmax, method,
                                                xtol, ftol, gtol, verbose, progress_bar,
                                                stages, **kwargs)
        else:
            raise ValueError(f"Unknown backend: {backend}. Choose from 'loky', 'threading', or 'sequential'.")

    def _fit_grouped_threading(self, data, params=None, wlmin: float = 1., wlmax: float = 6.,
                               method: str = "rietveld",
                               xtol: float = None, ftol: float = None, gtol: float = None,
                               verbose: bool = False,
                               progress_bar: bool = True,
                               stages: Optional[Union[str, Dict[str, Union[str, List[str]]]]] = None,
                               n_jobs: int = 10,
                               **kwargs):
        """
        Fit model to grouped data using threading backend (fallback when multiprocessing fails).

        Note: Threading doesn't provide true parallelism due to Python's GIL,
        but works when objects can't be pickled.
        """
        from joblib import Parallel, delayed
        import time

        try:
            from tqdm.auto import tqdm
        except ImportError:
            from tqdm import tqdm

        # Prepare fit arguments
        fit_kwargs = {
            'params': params,
            'wlmin': wlmin,
            'wlmax': wlmax,
            'method': method,
            'xtol': xtol,
            'ftol': ftol,
            'gtol': gtol,
            'verbose': verbose if verbose else False,
            'progress_bar': False,
            'stages': stages,
            **kwargs
        }

        def fit_single_group(idx):
            """Fit a single group using threading."""
            group_data = Data()
            group_data.table = data.groups[idx]
            group_data.L = data.L
            group_data.tstep = data.tstep

            try:
                result = self.fit(group_data, **fit_kwargs)
            except Exception as e:
                result = None
            return idx, result

        start_time = time.time()
        backend = 'threading'
        # print(f"Using threading backend (limited parallelism due to Python GIL)...")

        # Execute with threading
        if progress_bar:
            iterator = tqdm(data.indices, desc=f"Fitting {len(data.indices)} groups")
        else:
            iterator = data.indices

        results = Parallel(
            n_jobs=n_jobs,
            backend=backend,
            verbose=5 if verbose else 0
        )(delayed(fit_single_group)(idx) for idx in iterator)

        elapsed = time.time() - start_time
        if verbose:
            print(f"Completed in {elapsed:.2f}s using '{backend}' backend | {elapsed/len(data.indices):.3f}s per fit")

        # Collect results
        grouped_result = GroupedFitResult(group_shape=data.group_shape)
        failed_indices = []
        for idx, result in results:
            if result is not None:
                grouped_result.add_result(idx, result)
            else:
                failed_indices.append(idx)

        if failed_indices and verbose:
            warnings.warn(f"Fitting failed for {len(failed_indices)}/{len(data.indices)} groups. "
                         f"Failed indices: {failed_indices[:10]}{'...' if len(failed_indices) > 10 else ''}")

        return grouped_result

    def _fit_grouped_loky(self, data, params=None, wlmin: float = 1., wlmax: float = 6.,
                          method: str = "rietveld",
                          xtol: float = None, ftol: float = None, gtol: float = None,
                          verbose: bool = False,
                          progress_bar: bool = True,
                          stages: Optional[Union[str, Dict[str, Union[str, List[str]]]]] = None,
                          n_jobs: int = 10,
                          **kwargs):
        """
        Fit model to grouped data using true multiprocessing (ProcessPoolExecutor).

        This method provides true parallelism by:
        1. Serializing the model configuration to a pickleable dict
        2. Using ProcessPoolExecutor which reuses worker processes
        3. Each worker reconstructs the model and fits in parallel
        4. Returns only pickleable results

        Note: First batch has initialization overhead (~3s per worker for NCrystal),
        but ProcessPoolExecutor reuses workers so subsequent tasks are fast.
        """
        from concurrent.futures import ProcessPoolExecutor
        import time

        try:
            from tqdm.auto import tqdm
        except ImportError:
            from tqdm import tqdm

        # Serialize model configuration to a dict (no NCrystal objects)
        model_dict = self._to_dict()

        # Prepare fit arguments (must be pickleable - no lmfit.Parameters object)
        fit_kwargs = {
            'params': None,  # Will use model's params
            'wlmin': wlmin,
            'wlmax': wlmax,
            'method': method,
            'xtol': xtol,
            'ftol': ftol,
            'gtol': gtol,
            'verbose': False,  # Disable per-worker verbose
            'progress_bar': False,  # Disable per-worker progress bar
            'stages': stages,
        }
        # Add kwargs that are pickleable
        for k, v in kwargs.items():
            try:
                pickle.dumps(v)
                fit_kwargs[k] = v
            except (TypeError, pickle.PicklingError):
                if verbose:
                    print(f"Warning: kwarg '{k}' cannot be pickled, skipping")

        # Prepare worker arguments: (idx, model_dict, table_dict, L, tstep, fit_kwargs)
        worker_args = []
        for idx in data.indices:
            table_dict = data.groups[idx].to_dict()
            worker_args.append((idx, model_dict, table_dict, data.L, data.tstep, fit_kwargs))

        start_time = time.time()
        n_workers = min(n_jobs, len(data.indices))

        if progress_bar and verbose:
            print(f"Fitting {len(data.indices)} groups using multiprocessing (n_workers={n_workers})...")

        # Use ProcessPoolExecutor for true multiprocessing with worker reuse
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            if progress_bar:
                # Use tqdm for progress tracking
                results = list(tqdm(
                    executor.map(_fit_single_group_worker, worker_args),
                    total=len(worker_args),
                    desc="Fitting groups"
                ))
            else:
                results = list(executor.map(_fit_single_group_worker, worker_args))

        elapsed = time.time() - start_time
        if verbose:
            print(f"Completed in {elapsed:.2f}s using multiprocessing | {elapsed/len(data.indices):.3f}s per fit")

        # Collect results and reconstruct result objects
        grouped_result = GroupedFitResult(group_shape=data.group_shape)
        failed_indices = []
        error_messages = []
        for idx, result_dict in results:
            if result_dict is not None and 'error' not in result_dict:
                # Reconstruct result object from dict
                result = _reconstruct_result_from_dict(result_dict, model=self)
                grouped_result.add_result(idx, result)
            else:
                failed_indices.append(idx)
                if result_dict and 'error' in result_dict:
                    error_messages.append(f"{idx}: {result_dict['error']}")

        if failed_indices and verbose:
            error_details = ""
            if error_messages:
                error_details = "\n" + "\n".join(error_messages[:3])
                if len(error_messages) > 3:
                    error_details += f"\n... and {len(error_messages) - 3} more errors"
            warnings.warn(f"Fitting failed for {len(failed_indices)}/{len(data.indices)} groups. "
                         f"Failed indices: {failed_indices[:10]}{'...' if len(failed_indices) > 10 else ''}{error_details}")

        return grouped_result

    def _fit_grouped_sequential(self, data, params=None, wlmin: float = 1., wlmax: float = 6.,
                                 method: str = "rietveld",
                                 xtol: float = None, ftol: float = None, gtol: float = None,
                                 verbose: bool = False,
                                 progress_bar: bool = True,
                                 stages: Optional[Union[str, Dict[str, Union[str, List[str]]]]] = None,
                                 **kwargs):
        """
        Fit model to grouped data sequentially (no parallelization).

        This is useful for debugging or when parallel execution causes issues.
        """
        import time

        try:
            from tqdm.auto import tqdm
        except ImportError:
            from tqdm import tqdm

        # Prepare fit arguments
        fit_kwargs = {
            'params': params,
            'wlmin': wlmin,
            'wlmax': wlmax,
            'method': method,
            'xtol': xtol,
            'ftol': ftol,
            'gtol': gtol,
            'verbose': verbose,
            'progress_bar': False,
            'stages': stages,
            **kwargs
        }

        start_time = time.time()
        grouped_result = GroupedFitResult(group_shape=data.group_shape)
        failed_indices = []

        iterator = tqdm(data.indices, desc=f"Fitting {len(data.indices)} groups") if progress_bar else data.indices

        for idx in iterator:
            group_data = Data()
            group_data.table = data.groups[idx]
            group_data.L = data.L
            group_data.tstep = data.tstep

            try:
                result = self.fit(group_data, **fit_kwargs)
                grouped_result.add_result(idx, result)
            except Exception as e:
                failed_indices.append(idx)
                if verbose:
                    print(f"Fitting failed for group {idx}: {e}")

        elapsed = time.time() - start_time
        if verbose:
            print(f"Completed in {elapsed:.2f}s using 'sequential' backend | {elapsed/len(data.indices):.3f}s per fit")

        if failed_indices and verbose:
            warnings.warn(f"Fitting failed for {len(failed_indices)}/{len(data.indices)} groups. "
                         f"Failed indices: {failed_indices[:10]}{'...' if len(failed_indices) > 10 else ''}")

        return grouped_result

    def _create_stages_summary_table_enhanced(self, stage_results, resolved_param_groups, stage_names=None,
                                            method="rietveld", color=True):
        import pandas as pd

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
