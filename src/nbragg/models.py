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
import json
import os
import pickle


def _fit_single_group_worker(args):
    """
    Worker function for parallel fitting of a single group.

    This function is defined at module level so it can be pickled for multiprocessing.
    It reconstructs the model from a serialized dict, fits the data, and returns
    only pickleable results.

    Parameters
    ----------
    args : tuple
        (idx, model_dict, table_dict, L, tstep, fit_kwargs)

    Returns
    -------
    tuple
        (idx, result_dict) where result_dict contains pickleable fit results,
        or (idx, error_string) if fitting failed.
    """
    idx, model_dict, table_dict, L, tstep, fit_kwargs = args

    try:
        # Reconstruct model from dict (creates new NCrystal objects in this process)
        from nbragg.models import TransmissionModel
        from nbragg.data import Data
        import pandas as pd

        model = TransmissionModel._from_dict(model_dict)

        # Reconstruct data
        group_data = Data()
        group_data.table = pd.DataFrame(table_dict)
        group_data.L = L
        group_data.tstep = tstep

        # Fit
        result = model.fit(group_data, **fit_kwargs)

        # Extract only pickleable attributes
        result_dict = _extract_pickleable_result(result)

        return idx, result_dict

    except Exception as e:
        import traceback
        error_msg = f"{type(e).__name__}: {str(e)}\n{traceback.format_exc()}"
        return idx, {'error': error_msg}


def _extract_pickleable_result(fit_result):
    """
    Extract only pickleable attributes from a fit result.

    Returns a dictionary with all the important fit result data
    that can be safely passed between processes.
    """
    result_dict = {}

    # Core fit statistics
    for attr in ['success', 'chisqr', 'redchi', 'aic', 'bic',
                 'nvarys', 'ndata', 'nfev', 'message', 'method']:
        if hasattr(fit_result, attr):
            result_dict[attr] = getattr(fit_result, attr)

    # Parameters - serialize to JSON string
    if hasattr(fit_result, 'params'):
        result_dict['params_json'] = fit_result.params.dumps()

    # Best fit values
    if hasattr(fit_result, 'best_values'):
        result_dict['best_values'] = dict(fit_result.best_values)

    # Init values
    if hasattr(fit_result, 'init_values'):
        result_dict['init_values'] = dict(fit_result.init_values)

    # Residual and best_fit arrays
    if hasattr(fit_result, 'residual') and fit_result.residual is not None:
        result_dict['residual'] = fit_result.residual.tolist()
    if hasattr(fit_result, 'best_fit') and fit_result.best_fit is not None:
        result_dict['best_fit'] = fit_result.best_fit.tolist()

    # Covariance matrix
    if hasattr(fit_result, 'covar') and fit_result.covar is not None:
        result_dict['covar'] = fit_result.covar.tolist()

    # Variable names and init_vals
    if hasattr(fit_result, 'var_names'):
        result_dict['var_names'] = list(fit_result.var_names)
    if hasattr(fit_result, 'init_vals'):
        result_dict['init_vals'] = list(fit_result.init_vals)

    # User keywords (needed for plotting)
    if hasattr(fit_result, 'userkws') and fit_result.userkws:
        result_dict['userkws'] = dict(fit_result.userkws)

    # Data array (needed for plotting)
    if hasattr(fit_result, 'data') and fit_result.data is not None:
        result_dict['data'] = fit_result.data.tolist()

    # Weights array (needed for plotting)
    if hasattr(fit_result, 'weights') and fit_result.weights is not None:
        result_dict['weights'] = fit_result.weights.tolist()

    # Stage results if present (for rietveld/staged fits)
    if hasattr(fit_result, 'fit_stages') and fit_result.fit_stages:
        result_dict['fit_stages'] = [_extract_pickleable_result(stage) for stage in fit_result.fit_stages]

    return result_dict


def _reconstruct_result_from_dict(result_dict, model=None):
    """
    Reconstruct a minimal ModelResult-like object from a pickleable dict.

    This creates an object with the same interface as lmfit.ModelResult
    but from serialized data.
    """
    import lmfit
    import numpy as np

    # Create a minimal result object using SimpleNamespace
    from types import SimpleNamespace
    result = SimpleNamespace()

    # Restore basic attributes
    for attr in ['success', 'chisqr', 'redchi', 'aic', 'bic',
                 'nvarys', 'ndata', 'nfev', 'message', 'method']:
        if attr in result_dict:
            setattr(result, attr, result_dict[attr])

    # Restore parameters
    if 'params_json' in result_dict:
        result.params = lmfit.Parameters()
        result.params.loads(result_dict['params_json'])

    # Restore best values and init values
    if 'best_values' in result_dict:
        result.best_values = result_dict['best_values']
    if 'init_values' in result_dict:
        result.init_values = result_dict['init_values']

    # Restore arrays
    if 'residual' in result_dict:
        result.residual = np.array(result_dict['residual'])
    if 'best_fit' in result_dict:
        result.best_fit = np.array(result_dict['best_fit'])
    if 'covar' in result_dict:
        result.covar = np.array(result_dict['covar'])

    # Restore variable info
    if 'var_names' in result_dict:
        result.var_names = result_dict['var_names']
    if 'init_vals' in result_dict:
        result.init_vals = result_dict['init_vals']

    # Restore userkws (needed for plotting)
    if 'userkws' in result_dict:
        result.userkws = result_dict['userkws']

    # Restore data (needed for plotting)
    if 'data' in result_dict:
        result.data = np.array(result_dict['data'])

    # Restore weights (needed for plotting)
    if 'weights' in result_dict:
        result.weights = np.array(result_dict['weights'])

    # Restore stage results if present
    if 'fit_stages' in result_dict:
        result.fit_stages = [_reconstruct_result_from_dict(stage) for stage in result_dict['fit_stages']]

    # Add model reference if provided
    result.model = model

    return result


def _add_save_method_to_result(result):
    """
    Add a save() and fit_report() method to an lmfit.ModelResult object.

    This function monkey-patches the result object to add save functionality
    and a fit_report() method that returns the HTML representation.
    """
    def save(filename: str):
        """Save this fit result to a JSON file."""
        _save_result_impl(result, filename)

    def fit_report_html():
        """
        Return the HTML fit report for display in Jupyter notebooks.

        This method provides the same output as the automatic display
        when the result object is shown in a Jupyter cell.

        Returns:
        --------
        str
            HTML string containing the formatted fit results from lmfit.

        Examples:
        ---------
        >>> result = model.fit(data)
        >>> html_report = result.fit_report()
        >>> # Display in Jupyter:
        >>> from IPython.display import HTML, display
        >>> display(HTML(html_report))
        """
        if hasattr(result, '_repr_html_'):
            return result._repr_html_()
        else:
            # Fallback: return empty string if _repr_html_ is not available
            return ""

    # Store original fit_report if it exists
    original_fit_report = result.fit_report if hasattr(result, 'fit_report') else None

    # Add the methods to the result instance
    result.save = save
    result.fit_report = fit_report_html
    # Preserve access to original text-based fit_report
    if original_fit_report:
        result.fit_report_text = original_fit_report
    return result


def _save_result_impl(result, filename: str):
    """Implementation of result saving logic."""
    # Prepare fit result state
    state = {
        'version': '1.0',
        'class': 'ModelResult',
        'params': result.params.dumps(),
        'init_params': result.init_params.dumps() if hasattr(result, 'init_params') and result.init_params else None,
        'success': result.success if hasattr(result, 'success') else None,
        'message': result.message if hasattr(result, 'message') else None,
        'method': result.method if hasattr(result, 'method') else 'unknown',
        'nfev': result.nfev if hasattr(result, 'nfev') else None,
        'nvarys': result.nvarys if hasattr(result, 'nvarys') else None,
        'ndata': result.ndata if hasattr(result, 'ndata') else None,
        'nfree': result.nfree if hasattr(result, 'nfree') else None,
        'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
        'redchi': result.redchi if hasattr(result, 'redchi') else None,
        'aic': result.aic if hasattr(result, 'aic') else None,
        'bic': result.bic if hasattr(result, 'bic') else None,
    }

    # Save the fit result
    with open(filename, 'w') as f:
        json.dump(state, f, indent=2)

    # Save the model with fitted parameters
    model_filename = filename.replace('.json', '_model.json')
    if model_filename == filename:
        model_filename = filename.replace('.json', '') + '_model.json'

    if hasattr(result, 'model') and isinstance(result.model, TransmissionModel):
        # Temporarily update model params with fitted values
        original_params = result.model.params.copy()
        result.model.params = result.params

        # Save the model with fitted parameters
        result.model.save(model_filename)

        # Restore original params
        result.model.params = original_params


class GroupedFitResult:
    """
    Container for fit results from grouped data.

    Stores multiple ModelResult objects indexed by their group identifiers
    (integers, tuples, or strings depending on the data structure).

    Attributes:
    -----------
    results : dict
        Dictionary mapping group indices to lmfit.ModelResult objects.
    indices : list
        List of group indices in order.
    group_shape : tuple or None
        Shape of the grouped data ((ny, nx) for 2D, (n,) for 1D, None for named).

    Examples:
    ---------
    >>> # Access individual results
    >>> grouped_result = model.fit(grouped_data)
    >>> result_0_0 = grouped_result[(0, 0)]
    >>> result_0_0.plot()

    >>> # Plot parameter map
    >>> grouped_result.plot_parameter_map("thickness")

    >>> # Print summary
    >>> grouped_result.summary()
    """

    def __init__(self, group_shape=None):
        """
        Initialize an empty GroupedFitResult.

        Parameters:
        -----------
        group_shape : tuple or None
            Shape of the grouped data.
        """
        self.results = {}
        self.indices = []
        self.group_shape = group_shape

    def _normalize_index(self, index):
        """
        Normalize index for consistent lookup.
        Converts tuples to strings without spaces: (10, 20) -> "(10,20)"
        Accepts both "(10,20)" and "(10, 20)" string formats.
        """
        if isinstance(index, tuple):
            return str(index).replace(" ", "")
        elif isinstance(index, str):
            return index.replace(" ", "")
        else:
            return str(index)

    def _parse_string_index(self, string_idx):
        """
        Parse a string index back to its original form.
        "(10,20)" or "(10, 20)" -> (10, 20)
        "5" -> 5
        "center" -> "center"
        """
        import ast
        try:
            parsed = ast.literal_eval(string_idx)
            return parsed
        except (ValueError, SyntaxError):
            return string_idx

    def add_result(self, index, result):
        """
        Add a fit result for a specific group.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index.
        result : lmfit.ModelResult
            The fit result for this group.
        """
        # Normalize index for consistent storage
        normalized_index = self._normalize_index(index)
        self.results[normalized_index] = result
        if normalized_index not in self.indices:
            self.indices.append(normalized_index)

    def __getitem__(self, index):
        """
        Access a specific group's result.

        Supports flexible index access:
        - Tuples: (0, 0) or "(0,0)" or "(0, 0)"
        - Integers: 5 or "5"
        - Strings: "groupname"
        """
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise KeyError(f"Index {index} not found in results. Available: {self.indices}")
        return self.results[normalized_index]

    def __len__(self):
        """Return number of group results."""
        return len(self.results)

    def __repr__(self):
        """String representation."""
        return f"GroupedFitResult({len(self.results)} groups, shape={self.group_shape})"

    def plot(self, index, **kwargs):
        """
        Plot a specific group's fit result.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to plot.
            - For 2D grids: can use tuple (0, 0) or string "(0,0)" or "(0, 0)"
            - For 1D arrays: can use int 5 or string "5"
            - For named groups: use string "groupname"
        **kwargs
            Additional plotting parameters passed to result.plot().

        Returns:
        --------
        matplotlib.Axes
            The plot axes.
        """
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        # Get the individual fit result
        fit_result = self.results[normalized_index]

        # Check if this is a proper ModelResult or a SimpleNamespace (from loaded file)
        from types import SimpleNamespace
        if isinstance(fit_result, SimpleNamespace) or not hasattr(fit_result, 'userkws'):
            # Try to get the model from the result
            model = getattr(fit_result, 'model', None)

            # If no model available, raise an error
            if model is None:
                raise AttributeError(
                    f"Cannot plot index {index}: result was loaded from file without model information. "
                    "This typically happens with compact results. Try loading with the full model:\n"
                    "  model = TransmissionModel.load('model.json')\n"
                    "  result = GroupedFitResult.load('result.json', model=model)\n"
                    "  result.plot(index=...)"
                )

            # We have a model but missing userkws - this is an old saved file
            # Create empty userkws to allow plotting
            if not hasattr(fit_result, 'userkws'):
                fit_result.userkws = {}

        # Call the model's plot method but temporarily set fit_result to the correct one
        # This is needed because ModelResult.plot() delegates to Model.plot() where
        # self is the shared Model instance, not the individual ModelResult
        model = fit_result.model
        original_fit_result = getattr(model, 'fit_result', None)
        try:
            model.fit_result = fit_result
            return model.plot(**kwargs)
        finally:
            # Restore original fit_result
            if original_fit_result is not None:
                model.fit_result = original_fit_result
            elif hasattr(model, 'fit_result'):
                delattr(model, 'fit_result')

    def plot_parameter_map(self, param_name, query=None, kind=None, **kwargs):
        """
        Plot spatial map of a fitted parameter value, error, or fit statistic.

        Parameters:
        -----------
        param_name : str
            Name of the parameter to visualize. Can be:
            - Parameter name: "thickness", "norm", etc.
            - Parameter error: "thickness_err", "norm_err", etc.
            - Fit statistic: "redchi", "chisqr", "aic", "bic"
        query : str, optional
            Pandas query string to filter results (e.g., "redchi < 2").
            Can reference any parameter name, parameter_err, or statistic.
        kind : str, optional
            Plot type. If None (default), auto-detected based on group_shape:
            - 2D data: 'pcolormesh'
            - 1D data: 'line'
            - Named groups: 'bar'
            For 1D/named data, can also specify: 'line', 'bar', or 'errorbar'.
            Ignored for 2D data (always uses pcolormesh).
        **kwargs : dict, optional
            Additional plotting parameters:
            - cmap : str, optional
              Colormap for 2D maps (default: 'viridis').
            - title : str, optional
              Plot title (default: auto-generated).
            - vmin, vmax : float, optional
              Color scale limits.
            - figsize : tuple, optional
              Figure size (width, height) in inches.

        Returns:
        --------
        matplotlib.Axes
            The plot axes.

        Examples:
        ---------
        >>> # Plot parameter value
        >>> grouped_result.plot_parameter_map("thickness")
        >>>
        >>> # Plot parameter error
        >>> grouped_result.plot_parameter_map("thickness_err")
        >>>
        >>> # Plot redchi map
        >>> grouped_result.plot_parameter_map("redchi", cmap='hot')
        >>>
        >>> # Filter by redchi
        >>> grouped_result.plot_parameter_map("thickness", query="redchi < 2")
        >>>
        >>> # Complex query
        >>> grouped_result.plot_parameter_map("norm", query="redchi < 1.5 and thickness_err < 0.01")
        >>>
        >>> # 1D plot with error bars
        >>> grouped_result.plot_parameter_map("thickness", kind='errorbar')
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        # Auto-detect plot kind based on group_shape if not specified
        if kind is None:
            if self.group_shape and len(self.group_shape) == 2:
                kind = 'pcolormesh'
            elif self.group_shape and len(self.group_shape) == 1:
                kind = 'line'
            else:
                # Named groups or no shape
                kind = 'bar'

        # Build DataFrame with all parameters, errors, and statistics
        data_for_query = []
        param_values = {}

        for idx in self.indices:
            result = self.results[idx]
            row = {'index': idx}

            try:
                # Add all parameter values and errors
                for pname in result.params:
                    param = result.params[pname]
                    row[pname] = param.value
                    row[f"{pname}_err"] = param.stderr if param.stderr is not None else np.nan

                # Add fit statistics
                row['redchi'] = result.redchi if hasattr(result, 'redchi') else np.nan
                row['chisqr'] = result.chisqr if hasattr(result, 'chisqr') else np.nan
                row['aic'] = result.aic if hasattr(result, 'aic') else np.nan
                row['bic'] = result.bic if hasattr(result, 'bic') else np.nan
                row['nfev'] = result.nfev if hasattr(result, 'nfev') else np.nan

                data_for_query.append(row)

                # Extract the specific parameter value requested
                if param_name.endswith('_err'):
                    # Error requested
                    base_param = param_name[:-4]
                    if base_param in result.params:
                        param_values[idx] = result.params[base_param].stderr
                elif param_name in ['redchi', 'chisqr', 'aic', 'bic', 'nfev']:
                    # Statistic requested
                    param_values[idx] = getattr(result, param_name, np.nan)
                elif param_name in result.params:
                    # Parameter value requested
                    param_values[idx] = result.params[param_name].value
                else:
                    param_values[idx] = np.nan

            except Exception as e:
                param_values[idx] = np.nan

        # Apply query filter if provided
        indices_to_plot = self.indices
        if query:
            df = pd.DataFrame(data_for_query)
            try:
                filtered_df = df.query(query)
                indices_to_plot = [row['index'] for _, row in filtered_df.iterrows()]
                # Mask out filtered indices
                for idx in self.indices:
                    if idx not in indices_to_plot:
                        param_values[idx] = np.nan
            except Exception as e:
                print(f"Warning: Query failed: {e}")
                print("Plotting all data without filtering.")

        # Extract kwargs
        cmap = kwargs.pop("cmap", "viridis")
        title = kwargs.pop("title", None)
        vmin = kwargs.pop("vmin", None)
        vmax = kwargs.pop("vmax", None)
        figsize = kwargs.pop("figsize", None)

        # Create visualization based on group_shape
        if self.group_shape and len(self.group_shape) == 2:
            # 2D pcolormesh for proper block sizing
            # Extract unique x and y coordinates by parsing string indices
            xs = []
            ys = []
            for idx_str in self.indices:
                idx = self._parse_string_index(idx_str)
                if isinstance(idx, tuple) and len(idx) == 2:
                    xs.append(idx[0])
                    ys.append(idx[1])
            xs = sorted(set(xs))
            ys = sorted(set(ys))

            # Handle empty indices
            if len(xs) == 0 or len(ys) == 0:
                raise ValueError("No valid 2D indices found for plotting")

            # Calculate grid spacing (block size)
            x_spacing = xs[1] - xs[0] if len(xs) > 1 else 1
            y_spacing = ys[1] - ys[0] if len(ys) > 1 else 1

            # Create coordinate arrays including edges for pcolormesh
            x_edges = np.array(xs) - x_spacing / 2
            x_edges = np.append(x_edges, xs[-1] + x_spacing / 2)
            y_edges = np.array(ys) - y_spacing / 2
            y_edges = np.append(y_edges, ys[-1] + y_spacing / 2)

            # Create 2D array for values
            param_array = np.full((len(ys), len(xs)), np.nan)

            # Map indices to array positions
            x_map = {x: i for i, x in enumerate(xs)}
            y_map = {y: i for i, y in enumerate(ys)}

            for idx_str in self.indices:
                idx = self._parse_string_index(idx_str)
                if isinstance(idx, tuple) and len(idx) == 2:
                    x, y = idx
                    if x in x_map and y in y_map:
                        param_array[y_map[y], x_map[x]] = param_values[idx_str]

            fig, ax = plt.subplots(figsize=figsize)
            im = ax.pcolormesh(x_edges, y_edges, param_array, cmap=cmap, vmin=vmin, vmax=vmax,
                              shading='flat', **kwargs)
            ax.set_xlabel("X coordinate")
            ax.set_ylabel("Y coordinate")
            ax.set_aspect('equal')
            if title is None:
                title = f"{param_name} Map"
            ax.set_title(title)
            plt.colorbar(im, ax=ax, label=param_name)
            return ax

        elif self.group_shape and len(self.group_shape) == 1:
            # 1D plot with various styles - parse string indices back to ints
            indices_array = np.array([self._parse_string_index(idx) for idx in self.indices])
            # Replace None with np.nan for plotting
            values = np.array([param_values[idx] if param_values[idx] is not None else np.nan for idx in self.indices])

            # Get errors if available for errorbar plot
            errors = None
            if kind == 'errorbar' and not param_name.endswith('_err'):
                # Try to get errors for the parameter
                errors = []
                for idx in self.indices:
                    result = self.results[idx]
                    if param_name in result.params:
                        stderr = result.params[param_name].stderr
                        errors.append(stderr if stderr is not None else 0)
                    else:
                        errors.append(0)
                errors = np.array(errors)

            fig, ax = plt.subplots(figsize=figsize)

            if kind == 'line':
                ax.plot(indices_array, values, 'o-', **kwargs)
            elif kind == 'bar':
                ax.bar(indices_array, values, **kwargs)
            elif kind == 'errorbar':
                if errors is not None:
                    ax.errorbar(indices_array, values, yerr=errors, fmt='o-', capsize=5, **kwargs)
                else:
                    ax.plot(indices_array, values, 'o-', **kwargs)
            else:
                raise ValueError(f"Unknown kind '{kind}'. Must be 'line', 'bar', or 'errorbar'.")

            ax.set_xlabel("Index")
            ax.set_ylabel(param_name)
            if title is None:
                title = f"{param_name} vs Index"
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            return ax

        else:
            # Named indices - bar or line plot
            fig, ax = plt.subplots(figsize=figsize)
            positions = np.arange(len(self.indices))
            # Replace None with np.nan for plotting
            values = [param_values[idx] if param_values[idx] is not None else np.nan for idx in self.indices]

            # Get errors if available for errorbar plot
            errors = None
            if kind == 'errorbar' and not param_name.endswith('_err'):
                errors = []
                for idx in self.indices:
                    result = self.results[idx]
                    if param_name in result.params:
                        stderr = result.params[param_name].stderr
                        errors.append(stderr if stderr is not None else 0)
                    else:
                        errors.append(0)

            if kind == 'line':
                ax.plot(positions, values, 'o-', **kwargs)
            elif kind == 'bar':
                ax.bar(positions, values, **kwargs)
            elif kind == 'errorbar':
                if errors is not None:
                    ax.errorbar(positions, values, yerr=errors, fmt='o', capsize=5, **kwargs)
                else:
                    ax.plot(positions, values, 'o', **kwargs)
            else:
                raise ValueError(f"Unknown kind '{kind}'. Must be 'line', 'bar', or 'errorbar'.")

            ax.set_xticks(positions)
            ax.set_xticklabels(self.indices, rotation=45, ha='right')
            ax.set_ylabel(param_name)
            if title is None:
                title = f"{param_name} by Group"
            ax.set_title(title)
            plt.tight_layout()
            return ax

    def summary(self):
        """
        Print summary statistics for all group fits.

        Returns a pandas DataFrame with fit statistics and parameter values/errors for each group.

        Returns:
        --------
        pandas.DataFrame
            Summary table with columns: index, success, redchi, parameters, and parameter errors.
        """
        import pandas as pd
        import numpy as np

        summary_data = []
        for idx in self.indices:
            result = self.results[idx]
            row = {
                'index': str(idx),
                'success': result.success if hasattr(result, 'success') else None,
                'redchi': result.redchi if hasattr(result, 'redchi') else None,
                'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
                'nfev': result.nfev if hasattr(result, 'nfev') else None,
                'nvarys': result.nvarys if hasattr(result, 'nvarys') else None,
            }

            # Add all parameter values and errors
            if hasattr(result, 'params'):
                for param_name in result.params:
                    param = result.params[param_name]
                    row[param_name] = param.value
                    row[f"{param_name}_err"] = param.stderr if param.stderr is not None else np.nan

            summary_data.append(row)

        df = pd.DataFrame(summary_data)
        print("\nGrouped Fit Results Summary")
        print("=" * 80)
        print(df.to_string(index=False))
        print("=" * 80)
        return df

    def stages_summary(self, index):
        """
        Get the stages summary table for a specific group.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to get stages summary for.
            - For 2D grids: can use tuple (0, 0) or string "(0, 0)"
            - For 1D arrays: can use int 5 or string "5"
            - For named groups: use string "groupname"

        Returns:
        --------
        pandas.DataFrame or None
            The stages summary table for the specified group, or None if not available.

        Examples:
        ---------
        >>> result.stages_summary(index="empty2")
        >>> result.stages_summary(index=(0, 0))
        """
        # Normalize index for consistent lookup
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        result = self.results[normalized_index]

        if hasattr(result, 'stages_summary'):
            return result.stages_summary
        else:
            print(f"Warning: No stages_summary available for index {index}")
            return None

    def plot_total_xs(self, index, **kwargs):
        """
        Plot the total cross-section for a specific group.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to plot.
            - For 2D grids: can use tuple (0, 0) or string "(0, 0)"
            - For 1D arrays: can use int 5 or string "5"
            - For named groups: use string "groupname"
        **kwargs
            Additional plotting parameters passed to result.plot_total_xs().
            See ModelResult.plot_total_xs() for available options.

        Returns:
        --------
        matplotlib.Axes or tuple of matplotlib.Axes
            The plot axes (or tuple of axes if plot_residuals=True).

        Examples:
        ---------
        >>> result.plot_total_xs(index="empty1")
        >>> result.plot_total_xs(index=(0, 0), plot_bg=True, plot_dspace=True)
        >>> result.plot_total_xs(index=5, split_phases=True)
        """
        # Normalize index for consistent lookup
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        # Get the individual fit result
        fit_result = self.results[normalized_index]

        # Check if this is a proper ModelResult or a SimpleNamespace (from loaded file)
        from types import SimpleNamespace
        if isinstance(fit_result, SimpleNamespace) or not hasattr(fit_result, 'userkws'):
            # Try to get the model from the result
            model = getattr(fit_result, 'model', None)

            # If no model available, raise an error
            if model is None:
                raise AttributeError(
                    f"Cannot plot index {index}: result was loaded from file without model information. "
                    "This typically happens with compact results. Try loading with the full model:\n"
                    "  model = TransmissionModel.load('model.json')\n"
                    "  result = GroupedFitResult.load('result.json', model=model)\n"
                    "  result.plot_total_xs(index=...)"
                )

            # We have a model but missing userkws - this is an old saved file
            # Create empty userkws to allow plotting
            if not hasattr(fit_result, 'userkws'):
                fit_result.userkws = {}

        # Temporarily set model.fit_result to the correct one for this plot
        # (same issue as in plot() method - ModelResult delegates to shared Model)
        model = fit_result.model
        original_fit_result = getattr(model, 'fit_result', None)
        try:
            model.fit_result = fit_result
            return model.plot_total_xs(**kwargs)
        finally:
            # Restore original fit_result
            if original_fit_result is not None:
                model.fit_result = original_fit_result
            elif hasattr(model, 'fit_result'):
                delattr(model, 'fit_result')

    def _repr_html_(self):
        """
        HTML representation for Jupyter notebooks.

        Returns a formatted table summarizing all grouped fit results,
        including fit statistics and parameter values with errors.
        """
        import pandas as pd
        import numpy as np

        # Collect summary data
        summary_data = []
        for idx in self.indices:
            result = self.results[idx]
            row = {
                'index': str(idx),
                'success': result.success if hasattr(result, 'success') else None,
                'redchi': result.redchi if hasattr(result, 'redchi') else None,
                'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
                'nfev': result.nfev if hasattr(result, 'nfev') else None,
                'nvarys': result.nvarys if hasattr(result, 'nvarys') else None,
            }

            # Add all parameter values and errors
            if hasattr(result, 'params'):
                for param_name in result.params:
                    param = result.params[param_name]
                    row[param_name] = param.value
                    row[f"{param_name}_err"] = param.stderr if param.stderr is not None else np.nan

            summary_data.append(row)

        df = pd.DataFrame(summary_data)

        # Create HTML with styling
        html = f"""
        <div style="max-width: 100%; overflow-x: auto;">
            <h3>Grouped Fit Results Summary</h3>
            <p><b>Number of groups:</b> {len(self.indices)}</p>
            <p><b>Group shape:</b> {self.group_shape if self.group_shape else 'Named groups'}</p>
            {df.to_html(index=False, classes='dataframe', border=0, float_format=lambda x: f'{x:.4g}')}
        </div>
        """

        return html

    def summary(self):
        """
        Display the HTML summary table for all grouped fit results.

        In Jupyter notebooks, automatically displays the HTML table.
        Outside Jupyter, returns the HTML string.

        Returns:
        --------
        str or IPython.display.HTML
            HTML string (outside Jupyter) or displayed HTML (in Jupyter).

        Examples:
        ---------
        >>> result = model.fit(grouped_data)
        >>> result.summary()  # Auto-displays in Jupyter
        """
        html = self._repr_html_()

        # Try to detect if we're in a Jupyter environment
        try:
            from IPython.display import HTML, display
            from IPython import get_ipython
            if get_ipython() is not None:
                # We're in IPython/Jupyter - display the HTML
                display(HTML(html))
                return None
        except ImportError:
            pass

        # Not in Jupyter - return the HTML string
        return html

    def fit_report(self, index):
        """
        Display the HTML fit report for a specific group.

        In Jupyter notebooks, automatically displays the HTML report.
        Outside Jupyter, returns the HTML string.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to get the fit report for.
            - For 2D grids: can use tuple (0, 0) or string "(0, 0)"
            - For 1D arrays: can use int 5 or string "5"
            - For named groups: use string "groupname"

        Returns:
        --------
        str or IPython.display.HTML or None
            HTML string (outside Jupyter), displayed HTML (in Jupyter), or None.

        Examples:
        ---------
        >>> result = model.fit(grouped_data)
        >>> result.fit_report(index=(0, 0))  # Auto-displays in Jupyter
        """
        normalized_index = self._normalize_index(index)
        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        fit_result = self.results[normalized_index]

        # Check if it's a proper ModelResult or a SimpleNamespace
        if hasattr(fit_result, '_repr_html_'):
            html = fit_result._repr_html_()
        else:
            # If it's a SimpleNamespace (from loaded file), create a basic HTML report
            import pandas as pd

            html = '<div style="max-width: 900px;">\n'
            html += f'<h3>Fit Report for Index: {index}</h3>\n'

            # Parameters table
            if hasattr(fit_result, 'params'):
                html += '<h4>Parameters:</h4>\n'
                param_data = []
                for pname, param in fit_result.params.items():
                    if hasattr(param, 'value'):
                        param_data.append({
                            'Parameter': pname,
                            'Value': f"{param.value:.6g}",
                            'Std Error': f"{param.stderr:.6g}" if param.stderr else 'N/A',
                            'Vary': param.vary
                        })
                df = pd.DataFrame(param_data)
                html += df.to_html(index=False)

            # Fit statistics
            html += '<h4>Fit Statistics:</h4>\n'
            stats_data = {
                'Reduced χ²': getattr(fit_result, 'redchi', 'N/A'),
                'χ²': getattr(fit_result, 'chisqr', 'N/A'),
                'Data points': getattr(fit_result, 'ndata', 'N/A'),
                'Variables': getattr(fit_result, 'nvarys', 'N/A'),
                'Function evals': getattr(fit_result, 'nfev', 'N/A'),
                'Success': getattr(fit_result, 'success', 'N/A'),
            }
            stats_df = pd.DataFrame(list(stats_data.items()), columns=['Statistic', 'Value'])
            html += stats_df.to_html(index=False)

            html += '</div>'

        # Try to detect if we're in a Jupyter environment and auto-display
        try:
            from IPython.display import HTML, display
            from IPython import get_ipython
            if get_ipython() is not None:
                # We're in IPython/Jupyter - display the HTML
                display(HTML(html))
                return None
        except ImportError:
            pass

        # Not in Jupyter - return the HTML string
        return html

    def save(self, filename: str, compact: bool = True, model_filename: str = None):
        """
        Save grouped fit results to a single JSON file.

        Parameters:
        -----------
        filename : str
            Path to the output JSON file.
        compact : bool, optional
            If True (default), save only essential data (params, errors, redchi2) to save memory.
            If False, save full fit results.
        model_filename : str, optional
            Path to save the model configuration. Only used if compact=False.
            If None, model is saved to filename.replace('.json', '_model.json').
        """
        import json

        # Prepare grouped results structure
        grouped_state = {
            'version': '1.0',
            'class': 'GroupedFitResult',
            'group_shape': self.group_shape,
            'indices': [str(idx) for idx in self.indices],  # Convert to strings for JSON
            'results': {}
        }

        # Save each result
        for idx in self.indices:
            result = self.results[idx]
            idx_str = str(idx)

            if compact:
                # Save only essential data for map plotting
                result_data = {
                    'compact': True,
                    'params': {},
                    'redchi': result.redchi if hasattr(result, 'redchi') else None,
                    'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
                    'success': result.success if hasattr(result, 'success') else None,
                }
                # Extract parameter values and errors
                for param_name in result.params:
                    param = result.params[param_name]
                    result_data['params'][param_name] = {
                        'value': param.value,
                        'stderr': param.stderr,
                        'vary': param.vary,
                    }
            else:
                # Save full result
                result_data = {
                    'compact': False,
                    'params': result.params.dumps(),
                    'init_params': result.init_params.dumps() if hasattr(result, 'init_params') else None,
                    'success': result.success if hasattr(result, 'success') else None,
                    'message': result.message if hasattr(result, 'message') else None,
                    'nfev': result.nfev if hasattr(result, 'nfev') else None,
                    'nvarys': result.nvarys if hasattr(result, 'nvarys') else None,
                    'ndata': result.ndata if hasattr(result, 'ndata') else None,
                    'nfree': result.nfree if hasattr(result, 'nfree') else None,
                    'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
                    'redchi': result.redchi if hasattr(result, 'redchi') else None,
                    'aic': result.aic if hasattr(result, 'aic') else None,
                    'bic': result.bic if hasattr(result, 'bic') else None,
                }

                # Save plotting data (userkws only - data/weights excluded for performance)
                # userkws needs to be serializable, so convert any numpy arrays to lists
                if hasattr(result, 'userkws') and result.userkws:
                    import numpy as np
                    serializable_userkws = {}
                    for key, value in result.userkws.items():
                        if isinstance(value, np.ndarray):
                            serializable_userkws[key] = value.tolist()
                        elif hasattr(value, '__iter__') and not isinstance(value, (str, dict)):
                            # Try to convert other iterables
                            try:
                                serializable_userkws[key] = list(value)
                            except:
                                # Skip non-serializable values
                                pass
                        else:
                            serializable_userkws[key] = value
                    result_data['userkws'] = serializable_userkws

                # Save stages_summary if available (for multi-stage fits)
                if hasattr(result, 'stages_summary') and result.stages_summary is not None:
                    import pandas as pd
                    if isinstance(result.stages_summary, pd.DataFrame):
                        result_data['stages_summary'] = result.stages_summary.to_json(orient='split')

            grouped_state['results'][idx_str] = result_data

        # Save to file
        with open(filename, 'w') as f:
            json.dump(grouped_state, f, indent=2)

        # Save model if not compact
        if not compact and model_filename != '' and len(self.indices) > 0:
            if model_filename is None:
                model_filename = filename.replace('.json', '_model.json')
                if model_filename == filename:
                    model_filename = filename.replace('.json', '') + '_model.json'

            # Get model from first result
            first_result = self.results[self.indices[0]]
            if hasattr(first_result, 'model') and isinstance(first_result.model, TransmissionModel):
                first_result.model.save(model_filename)

    @classmethod
    def load(cls, filename: str, model_filename: str = None, model: 'TransmissionModel' = None):
        """
        Load grouped fit results from a JSON file.

        Parameters:
        -----------
        filename : str
            Path to the saved JSON file.
        model_filename : str, optional
            Path to the model configuration file. Only needed if full results were saved.
        model : TransmissionModel, optional
            Existing model to use instead of loading from file.

        Returns:
        --------
        GroupedFitResult
            Loaded grouped fit results.
        """
        import json
        import ast

        with open(filename, 'r') as f:
            grouped_state = json.load(f)

        # Create new instance
        group_shape = tuple(grouped_state['group_shape']) if grouped_state['group_shape'] else None
        grouped_result = cls(group_shape=group_shape)

        # Parse indices back to original types
        indices_str = grouped_state['indices']
        indices = []
        indices_str_map = {}  # Map to keep track of original string representations
        for idx_str in indices_str:
            try:
                # Try to evaluate as tuple/int
                idx = ast.literal_eval(idx_str)
            except (ValueError, SyntaxError):
                # Keep as string
                idx = idx_str
            indices.append(idx)
            indices_str_map[idx if not isinstance(idx, (tuple, int)) else idx] = idx_str

        # Try to load model for compact results (for plotting support)
        model_for_compact = None
        if any(grouped_state['results'][idx_str]['compact'] for idx_str in indices_str):
            # At least one compact result - try to load model
            try:
                if model is None and model_filename is None:
                    model_filename = filename.replace('.json', '_model.json')
                    if model_filename == filename:
                        model_filename = filename.replace('.json', '') + '_model.json'
                if model is None:
                    model_for_compact = TransmissionModel.load(model_filename)
                else:
                    model_for_compact = model
            except (FileNotFoundError, Exception):
                # Model not available - compact results won't support plotting
                pass

        # Load each result
        for i, idx in enumerate(indices):
            idx_str = indices_str[i]  # Use original string representation
            result_data = grouped_state['results'][idx_str]

            if result_data['compact']:
                # Create a minimal result object for compact storage
                # This is a simple container, not a full ModelResult
                class CompactResult:
                    def __init__(self, data, model=None):
                        from lmfit import Parameters
                        self.params = Parameters()
                        for param_name, param_data in data['params'].items():
                            self.params.add(param_name,
                                          value=param_data['value'],
                                          vary=param_data['vary'])
                            self.params[param_name].stderr = param_data['stderr']
                        self.redchi = data['redchi']
                        self.chisqr = data['chisqr']
                        self.success = data['success']
                        self.compact = True
                        self.model = model  # Store reference to model for plotting

                result = CompactResult(result_data, model=model_for_compact)
            else:
                # Reconstruct full ModelResult
                from lmfit import Parameters, minimize
                from lmfit.model import ModelResult

                # Load or use provided model
                if model is None:
                    if model_filename is None:
                        model_filename = filename.replace('.json', '_model.json')
                        if model_filename == filename:
                            model_filename = filename.replace('.json', '') + '_model.json'
                    model = TransmissionModel.load(model_filename)

                # Create minimal result object
                params = Parameters()
                params.loads(result_data['params'])

                result = ModelResult(model, params)
                result.success = result_data['success']
                result.message = result_data['message']
                result.nfev = result_data['nfev']
                result.nvarys = result_data['nvarys']
                result.ndata = result_data['ndata']
                result.nfree = result_data['nfree']
                result.chisqr = result_data['chisqr']
                result.redchi = result_data['redchi']
                result.aic = result_data['aic']
                result.bic = result_data['bic']

                if result_data['init_params']:
                    init_params = Parameters()
                    init_params.loads(result_data['init_params'])
                    result.init_params = init_params

                # Restore userkws (data/weights not saved for performance)
                if 'userkws' in result_data:
                    result.userkws = result_data['userkws']
                else:
                    # Ensure userkws exists even if not in file (for old saves)
                    result.userkws = {}

                # Restore stages_summary if available (for multi-stage fits)
                if 'stages_summary' in result_data and result_data['stages_summary'] is not None:
                    import pandas as pd
                    from io import StringIO
                    result.stages_summary = pd.read_json(StringIO(result_data['stages_summary']), orient='split')

            grouped_result.add_result(idx, result)

        return grouped_result


class TransmissionModel(lmfit.Model):
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
            from nbragg.data import Data
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

    def _make_basic_params(self, vary=True):
        params = lmfit.Parameters()
        params.add("thickness", value=1., min=0., vary=vary)
        params.add("norm", value=1., min=0., vary=vary)
        return params

    def _make_temperature_params(self):
        params = lmfit.Parameters()
        params.add("temp", value=293.15, min=0., vary=False)  # Default temperature in Kelvin, always fixed by default
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
                    split_phases: bool = False,
                    plot_residuals: bool = False,
                    weight_label_position: str = 'right',
                    figsize: tuple = None,
                    height_ratios: list = None,
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
        split_phases: bool, optional
            If True, plots individual phase contributions with different colors and weight labels.
        plot_residuals: bool, optional
            If True, creates a 2-panel plot with residuals in the bottom panel.
        weight_label_position : str, optional
            Position of weight labels when split_phases=True. Options:
            - 'right': Labels on right edge of plot (default)
            - 'left': Labels above each curve on the left (respects log scale)
            - 'legend': Include weights in legend labels
            - 'none' or None: No weight labels
        figsize : tuple, optional
            Figure size as (width, height). Default is (6, 4) for single panel,
            (6, 5) for residuals panel.
        height_ratios : list, optional
            Height ratios for panels when plot_residuals=True.
            Default is [6, 1] (main panel 6x larger than residuals).
            Example: [3, 1] for different ratio.
        color : str, optional
            Color for the total cross section line. Default is "0.1" (dark gray).
        title : str, optional
            Plot title. Default is "Total Cross-Section: {material_name}".
        logy : bool, optional
            If True, use logarithmic scale for y-axis. Default is True.
        ylim : tuple, optional
            Y-axis limits as (ymin, ymax). Example: ylim=(1e-2, 20).
        xlim : tuple, optional
            X-axis limits as (xmin, xmax). Example: xlim=(1.0, 5.0).
        legend_loc : str, optional
            Legend location. Default is 'lower right' for data, 'best' otherwise.
        legend_fontsize : int, optional
            Legend font size. Default is 9.
        residuals_ylim : tuple, optional
            Y-axis limits for residuals panel. Default is (-2.5, 2.5).
        **kwargs : dict, optional
            Additional matplotlib parameters passed to ax.set() (e.g., xlabel, ylabel, etc.).

        Returns
        -------
        matplotlib.axes.Axes or tuple of Axes
            The axes of the plot. Returns tuple (ax_main, ax_residual) if plot_residuals=True.

        Notes
        -----
        This function generates a plot showing the total cross-section data and the best-fit curve.
        If `plot_bg` is True, it will also plot the background function.
        If `split_phases` is True, shows individual phase contributions with weights.
        Can be used both after fitting (using fit_result) or before fitting (using model params).
        """
        import matplotlib.pyplot as plt
        import numpy as np

        # Set default figsize and height_ratios
        if figsize is None:
            figsize = (6, 5) if plot_residuals else (6, 4)
        if height_ratios is None:
            height_ratios = [6, 1]

        # Create figure with or without residuals panel
        if plot_residuals:
            fig, (ax, ax_res) = plt.subplots(2, 1, figsize=figsize,
                                             height_ratios=height_ratios, sharex=True)
        else:
            fig, ax = plt.subplots(figsize=figsize)
            ax_res = None    

        # Determine which results to use
        data_xs = None
        data_err = None
        has_data = False

        if stage is not None and hasattr(self, "fit_stages") and self.fit_stages:
            # Use specific stage results
            if stage < 1 or stage > len(self.fit_stages):
                raise ValueError(f"Stage {stage} not available. Available stages: 1-{len(self.fit_stages)}")

            # Get stage results
            stage_result = self.fit_stages[stage - 1]  # Convert to 0-indexed
            params = stage_result.params
            wavelength = np.linspace(1.0, 10.0, 1000)  # Adjust range as needed
            fit_label = f"Stage {stage} fit"

        elif hasattr(self, "fit_result") and self.fit_result is not None:
            # Use final fit results
            wavelength = self.fit_result.userkws["wl"]
            params = self.fit_result.params
            fit_label = "Best fit"

            # Extract data for plotting
            if hasattr(self.fit_result, 'data') and hasattr(self.fit_result, 'weights'):
                trans_data = self.fit_result.data
                weights = self.fit_result.weights
                err_trans = 1.0 / weights

                # Convert transmission back to cross section
                # σ = -ln(T) / (thickness * n)
                thickness = params['thickness'].value
                norm = params.get('norm', params.get('normalization', type('obj', (), {'value': 1.0}))).value

                # Get background if it exists
                bg = np.zeros_like(trans_data)
                if self.background is not None:
                    bg = self.background.function(wl=wavelength, **params.valuesdict())

                # Calculate cross section from transmission: T = norm * exp(-σ * thickness * n) * (1 - bg) + bg
                # Rearranging: σ = -ln((T - bg) / (norm * (1 - bg))) / (thickness * n)
                trans_corrected = (trans_data - bg) / (norm * (1 - bg) + 1e-10)
                trans_corrected = np.clip(trans_corrected, 1e-10, 1.0)  # Ensure valid range

                # Get atomic density from model (same as used in transmission calculation)
                n = self.atomic_density

                data_xs = -np.log(trans_corrected) / (thickness * n)
                # Error propagation: Δσ = |dσ/dT| * ΔT = σ/T * ΔT (approximately)
                data_err = np.abs(data_xs * err_trans / (trans_data + 1e-10))
                has_data = True

        else:
            # Use model (no fit yet)
            fit_label = "Model"
            params = self.params
            wavelength = np.linspace(1.0, 10.0, 1000)  # Adjust range as needed

        # Calculate total cross section and individual phase contributions
        xs_total = self.cross_section(wavelength, **params.valuesdict())

        # Get individual phase cross sections if split_phases is True
        if split_phases:
            phase_xs = {}
            phase_weights = {}
            for phase in self.cross_section.phases:
                # Get cross section for this phase
                phase_xs[phase] = self.cross_section.get_phase_xs(wavelength, phase, **params.valuesdict())
                # Get weight
                if phase in params:
                    phase_weights[phase] = params[phase].value
                else:
                    # Try to find weight parameter
                    weight_found = False
                    for key in params:
                        if phase in key and 'weight' not in key:
                            phase_weights[phase] = params[key].value
                            weight_found = True
                            break
                    if not weight_found:
                        phase_weights[phase] = self.cross_section.materials.get(phase, {}).get('weight', 1.0)

        # Plot settings - extract specific parameters
        color = kwargs.pop("color", "0.1")
        title = kwargs.pop("title", f"Total Cross-Section: {self.cross_section.name}")
        logy = kwargs.pop("logy", True)  # Default to log scale
        ylim = kwargs.pop("ylim", None)
        xlim = kwargs.pop("xlim", None)
        legend_loc = kwargs.pop("legend_loc", None)
        legend_fontsize = kwargs.pop("legend_fontsize", 9)

        # Plot data if available
        if has_data and data_xs is not None:
            ax.errorbar(wavelength, data_xs, yerr=data_err,
                       marker=".", ms=2, ls="none",
                       color="#B0A1BA", ecolor="0.8",
                       zorder=-1, label="Data")

        # Plot individual phases if requested
        if split_phases:
            # Get colormap
            cmap = plt.cm.turbo
            n_phases = len(self.cross_section.phases)
            colors = cmap(np.linspace(0., 1, n_phases))

            # Sort phases by weight
            sorted_phases = sorted(phase_weights.items(), key=lambda x: x[1], reverse=True)
            total_weight = sum(phase_weights.values())

            # Plot each phase with appropriate labels
            for i, (phase, weight) in enumerate(sorted_phases):
                if phase in phase_xs:
                    weighted_xs = phase_xs[phase] * weight
                    percentage = (weight / total_weight) * 100

                    # Create label based on weight_label_position
                    if weight_label_position == 'legend':
                        phase_label = f"{phase}: {percentage:>3.1f}%"
                    else:
                        phase_label = f"{phase}"

                    ax.plot(wavelength, weighted_xs, lw=1, color=colors[i],
                           label=phase_label, zorder=5)

            # Add weight labels as text annotations (not in legend)
            if weight_label_position == 'right':
                # Labels on right edge of plot
                xlim_current = ax.get_xlim()
                x_label = xlim_current[1] * 0.98

                for i, (phase, weight) in enumerate(sorted_phases):
                    if phase in phase_xs:
                        weighted_xs = phase_xs[phase] * weight
                        y_pos = weighted_xs[-1]
                        # Filter out very small contributions
                        if y_pos > xs_total[-1] * 0.01:
                            percentage = (weight / total_weight) * 100
                            ax.text(x_label, y_pos, f"{phase}: {percentage:>3.1f}%",
                                   color=colors[i], fontsize=8, rotation=4,
                                   va='center', ha='right')

            elif weight_label_position == 'left':
                # Labels above each curve on the left
                xlim_current = ax.get_xlim()
                x_label = xlim_current[0] + (xlim_current[1] - xlim_current[0]) * 0.05

                for i, (phase, weight) in enumerate(sorted_phases):
                    if phase in phase_xs:
                        weighted_xs = phase_xs[phase] * weight
                        # Find y-position at the label x position
                        idx = np.argmin(np.abs(wavelength - x_label))
                        y_pos = weighted_xs[idx]
                        # Filter out very small contributions
                        if y_pos > xs_total[idx] * 0.01:
                            percentage = (weight / total_weight) * 100
                            # Offset label upward to avoid overlap with curve
                            # Use multiplicative offset for log scale, small additive for linear
                            if logy:
                                y_label = y_pos * 1.2  # 20% higher in log space
                            else:
                                # Use small multiplicative factor for linear scale too
                                y_label = y_pos * 1.05  # 5% higher

                            ax.text(x_label, y_label, f"{phase}: {percentage:>3.1f}%",
                                   color=colors[i], fontsize=8,
                                   va='bottom', ha='left')

        # Plot total cross-section
        ax.plot(wavelength, xs_total, color=color, label=fit_label, zorder=10, lw=1.5)
        ax.set_ylabel("Cross-Section [barn]")
        if not plot_residuals:
            ax.set_xlabel("λ [Å]")
        ax.set_title(title)

        # Apply log scale if requested
        if logy:
            ax.set_yscale("log")

        # Apply axis limits if provided
        if ylim is not None:
            ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)

        # Plot background if requested
        if plot_bg and self.background:
            bg = self.background.function(wl=wavelength, **params.valuesdict())
            ax.plot(wavelength, bg, color="orange", linestyle="--", label="Background")

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

        # Add legend
        if has_data:
            legend_title = f"$\\chi^2$: {self.fit_result.redchi:.2f}" if hasattr(self, 'fit_result') else None
            loc = legend_loc if legend_loc is not None else 'lower right'
            ax.legend(fontsize=legend_fontsize, loc=loc, title=legend_title,
                     title_fontsize=legend_fontsize, reverse=True)
        else:
            loc = legend_loc if legend_loc is not None else 'best'
            ax.legend(fontsize=legend_fontsize, loc=loc)

        # Plot residuals if requested
        if plot_residuals and has_data and data_xs is not None:
            residuals = data_xs - xs_total
            ax_res.plot(wavelength, residuals, marker=".", ms=1,
                       color="#B0A1BA", ls="none", zorder=-1)
            ax_res.axhline(0, color="0.1", lw=1, ls="-")
            ax_res.set_ylabel("Residuals [barn]", labelpad=13)
            ax_res.set_xlabel("λ [Å]")

            # Allow user to override residuals ylim
            residuals_ylim = kwargs.pop("residuals_ylim", [-2.5, 2.5])
            ax_res.set_ylim(residuals_ylim)
            plt.subplots_adjust(hspace=0.05)

        # Apply any remaining kwargs to the main axes
        # This allows users to pass additional matplotlib parameters
        if kwargs:
            try:
                ax.set(**kwargs)
            except:
                # If set() fails, ignore silently (kwargs might be for other purposes)
                pass

        plt.tight_layout()

        if plot_residuals:
            return ax, ax_res
        else:
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

    def save(self, filename: str):
        """
        Save the TransmissionModel configuration to a JSON file.

        This method saves all the necessary information to reconstruct the model,
        including cross-section materials, parameters, response and background types,
        and other configuration settings. It avoids pickling NCrystal objects by
        storing only the material specifications.

        Parameters
        ----------
        filename : str
            Path to the output JSON file.

        Notes
        -----
        The saved file can be loaded using TransmissionModel.load() to reconstruct
        the model with the same configuration.

        Examples
        --------
        >>> model = TransmissionModel(cross_section, vary_background=True)
        >>> model.save('my_model.json')
        """
        import json

        # Prepare model state dictionary
        state = {
            'version': '1.0',
            'class': 'TransmissionModel',
            'materials': self._materials,
            'cross_section_name': self.cross_section.name,
            'cross_section_total_weight': self.cross_section.total_weight,
            'cross_section_extinction': self.cross_section.extinction,
            'tof_length': self.tof_length,
            'params': self.params.dumps(),  # lmfit Parameters to JSON
            'stages': self._stages,
        }

        # Save response configuration if it exists
        if self.response is not None:
            state['response'] = {
                'kind': self.response.kind,
                'params': self.response.params.dumps()
            }
        else:
            state['response'] = None

        # Save background configuration if it exists
        if self.background is not None:
            # Infer background kind from parameters since Background doesn't store it
            bg_kind = 'polynomial3'  # default
            bg_params = list(self.background.params.keys())
            if 'k' in bg_params:
                bg_kind = 'sample_dependent'
            elif len([p for p in bg_params if p.startswith('bg')]) >= 5:
                bg_kind = 'polynomial5'
            elif len([p for p in bg_params if p.startswith('bg')]) == 3:
                bg_kind = 'polynomial3'
            elif len([p for p in bg_params if p.startswith('bg')]) == 1:
                bg_kind = 'constant'
            elif len(bg_params) == 0:
                bg_kind = 'none'

            state['background'] = {
                'kind': bg_kind,
                'params': self.background.params.dumps()
            }
        else:
            state['background'] = None

        # Write to file
        with open(filename, 'w') as f:
            json.dump(state, f, indent=2)

    def _to_dict(self) -> dict:
        """
        Convert model configuration to a pickleable dictionary.

        This is a fast alternative to save() that avoids file I/O.
        The resulting dict can be passed between processes and used
        to reconstruct the model with _from_dict().

        Returns
        -------
        dict
            A dictionary containing all model configuration needed for reconstruction.
            All values are pure Python types (no NCrystal objects).
        """
        # Process materials to use original specifications instead of virtual .nbragg files
        # This ensures the materials can be reconstructed in a subprocess
        materials_for_serialization = {}
        for name, mat_spec in self._materials.items():
            mat_copy = dict(mat_spec)
            # Use original material spec if available, otherwise keep current
            if '_original_mat' in mat_copy:
                mat_copy['mat'] = mat_copy['_original_mat']
            materials_for_serialization[name] = mat_copy

        # Prepare model state dictionary (same as save, but without file I/O)
        state = {
            'version': '1.0',
            'class': 'TransmissionModel',
            'materials': materials_for_serialization,
            'cross_section_name': self.cross_section.name,
            'cross_section_total_weight': self.cross_section.total_weight,
            'cross_section_extinction': self.cross_section.extinction,
            'tof_length': self.tof_length,
            'params': self.params.dumps(),  # lmfit Parameters to JSON string
            'stages': self._stages,
        }

        # Save response configuration if it exists
        if self.response is not None:
            state['response'] = {
                'kind': self.response.kind,
                'params': self.response.params.dumps()
            }
        else:
            state['response'] = None

        # Save background configuration if it exists
        if self.background is not None:
            # Infer background kind from parameters
            bg_kind = 'polynomial3'  # default
            bg_params = list(self.background.params.keys())
            if 'k' in bg_params:
                bg_kind = 'sample_dependent'
            elif len([p for p in bg_params if p.startswith('bg')]) >= 5:
                bg_kind = 'polynomial5'
            elif len([p for p in bg_params if p.startswith('bg')]) == 3:
                bg_kind = 'polynomial3'
            elif len([p for p in bg_params if p.startswith('bg')]) == 1:
                bg_kind = 'constant'
            elif len(bg_params) == 0:
                bg_kind = 'none'

            state['background'] = {
                'kind': bg_kind,
                'params': self.background.params.dumps()
            }
        else:
            state['background'] = None

        return state

    @classmethod
    def _from_dict(cls, state: dict) -> 'TransmissionModel':
        """
        Reconstruct a TransmissionModel from a dictionary.

        This is a fast alternative to load() that avoids file I/O.
        Used for parallel fitting where each worker needs to reconstruct
        the model from a pickleable configuration.

        Parameters
        ----------
        state : dict
            Dictionary from _to_dict() containing model configuration.

        Returns
        -------
        TransmissionModel
            A new TransmissionModel instance with the same configuration.
        """
        return cls._load_from_model(state)

    @classmethod
    def load(cls, filename: str):
        """
        Load a TransmissionModel from a JSON file (model or result).

        This method can load both model configuration files and fit result files.
        It automatically detects the file type and loads accordingly.

        Parameters
        ----------
        filename : str
            Path to the input JSON file (model or result).

        Returns
        -------
        TransmissionModel
            The reconstructed model with all parameters and settings restored.
            If loading from a result file, the model will have a `.result` attribute
            containing the loaded fit result.

        Notes
        -----
        The model is reconstructed by creating a new CrossSection from the saved
        material specifications and then initializing a new TransmissionModel with
        the saved parameters.

        When loading from a result file, the model is initialized with the fitted
        parameters and the result object is attached to `model.result`.

        Examples
        --------
        >>> # Load from model file
        >>> model = TransmissionModel.load('my_model.json')
        >>> result = model.fit(data)
        >>>
        >>> # Load from result file
        >>> model = TransmissionModel.load('my_result.json')
        >>> model.result.plot()  # Access the loaded result
        >>> print(model.result.redchi)
        """
        with open(filename, 'r') as f:
            state = json.load(f)

        # Verify version
        if state.get('version') != '1.0':
            warnings.warn(f"Loading file saved with version {state.get('version')}, "
                         f"current version is 1.0. Compatibility issues may occur.")

        # Detect file type
        if state.get('class') == 'ModelResult':
            # This is a result file
            return cls._load_from_result(filename, state)
        elif state.get('class') == 'TransmissionModel':
            # This is a model file
            return cls._load_from_model(state)
        else:
            raise ValueError(f"Unknown file type: {state.get('class')}")

    @classmethod
    def _load_from_model(cls, state):
        """Load a TransmissionModel from a model state dict."""
        # Reconstruct CrossSection
        cross_section = CrossSection(
            materials=state['materials'],
            name=state['cross_section_name'],
            total_weight=state['cross_section_total_weight']
        )

        # Restore extinction if it exists
        if 'cross_section_extinction' in state and state['cross_section_extinction']:
            cross_section.extinction = state['cross_section_extinction']

        # Load parameters
        params = lmfit.Parameters()
        params.loads(state['params'])

        # Determine response and background types
        response_kind = state['response']['kind'] if state['response'] is not None else None
        background_kind = state['background']['kind'] if state['background'] is not None else None

        # Create new model WITHOUT vary flags to avoid overwriting loaded params
        model = cls(
            cross_section=cross_section,
            params=params,
            response=response_kind if response_kind else "jorgensen",
            background=background_kind if background_kind else "polynomial3",
            tof_length=state['tof_length']
        )

        # Manually create response and background objects if they existed
        if response_kind is not None:
            model.response = Response(kind=response_kind, vary=False)
            # Update with loaded params
            for param_name in model.response.params.keys():
                if param_name in params:
                    model.response.params[param_name] = params[param_name]

        if background_kind is not None:
            model.background = Background(kind=background_kind, vary=False)
            # Store the kind attribute for consistency (Background class doesn't store it by default)
            model.background.kind = background_kind
            # Update with loaded params
            for param_name in model.background.params.keys():
                if param_name in params:
                    model.background.params[param_name] = params[param_name]

        # Restore stages
        model._stages = state['stages']

        return model

    @classmethod
    def _load_from_result(cls, filename, result_state):
        """Load a TransmissionModel from a result state dict and reconstruct the result."""
        # Load the associated model file
        model_filename = filename.replace('.json', '_model.json')
        if model_filename == filename:
            model_filename = filename.replace('.json', '') + '_model.json'

        if not os.path.exists(model_filename):
            raise FileNotFoundError(
                f"Model file {model_filename} not found. "
                f"Result files require an associated model file."
            )

        # Load the model
        with open(model_filename, 'r') as f:
            model_state = json.load(f)

        model = cls._load_from_model(model_state)

        # Reconstruct the result object
        result = cls._reconstruct_result(model, result_state)

        # Attach the result to the model
        model.result = result

        return model

    @classmethod
    def _reconstruct_result(cls, model, result_state):
        """
        Reconstruct a ModelResult object from saved state.

        This creates a "mock" ModelResult that has all the essential attributes
        and methods, including plot, _html_repr_, etc.
        """
        # Create a minimal result-like object
        result = lmfit.minimizer.MinimizerResult()

        # Load parameters
        result.params = lmfit.Parameters()
        result.params.loads(result_state['params'])

        if result_state['init_params'] is not None:
            result.init_params = lmfit.Parameters()
            result.init_params.loads(result_state['init_params'])
        else:
            result.init_params = None

        # Restore fit statistics
        result.success = result_state.get('success')
        result.message = result_state.get('message')
        result.nfev = result_state.get('nfev')
        result.nvarys = result_state.get('nvarys')
        result.ndata = result_state.get('ndata')
        result.nfree = result_state.get('nfree')
        result.chisqr = result_state.get('chisqr')
        result.redchi = result_state.get('redchi')
        result.aic = result_state.get('aic')
        result.bic = result_state.get('bic')

        # Add additional attributes that lmfit expects for _repr_html_
        result.method = result_state.get('method', 'loaded')
        result.aborted = False
        result.errorbars = True
        result.var_names = [name for name in result.params.keys() if result.params[name].vary]
        result.covar = None
        result.init_vals = result.init_params.valuesdict() if result.init_params else {}

        # Attach the model
        result.model = model

        # Add model-specific methods
        result.plot = model.plot
        result.plot_total_xs = model.plot_total_xs
        result.show_available_params = model.show_available_params

        if model.response is not None:
            result.response = model.response
            result.response.params = result.params

        if model.background is not None:
            result.background = model.background

        if hasattr(model, 'stages_summary'):
            result.stages_summary = model.stages_summary

        # Add the save method
        result = _add_save_method_to_result(result)

        return result


# Module-level functions for saving and loading fit results
def save_result(result, filename: str, model_filename: str = None):
    """
    Save a ModelResult (fit result) to JSON file(s).

    This function saves the fit results, including fitted parameters, statistics,
    and optionally the model configuration. It avoids the ctypes pickle issue by
    storing only serializable data.

    Parameters
    ----------
    result : lmfit.ModelResult
        The fit result object to save.
    filename : str
        Path to the output JSON file for the fit results.
    model_filename : str, optional
        Path to save the model configuration. If None, model is saved to
        filename.replace('.json', '_model.json'). If you don't want to save
        the model separately, pass an empty string ''.

    Notes
    -----
    The fit result can be loaded using load_result() to reconstruct both the
    model and the fit results.

    Examples
    --------
    >>> result = model.fit(data)
    >>> save_result(result, 'my_fit.json')
    >>> # Later...
    >>> loaded_result = load_result('my_fit.json')
    """
    import json

    # Prepare fit result state
    state = {
        'version': '1.0',
        'class': 'ModelResult',
        'params': result.params.dumps(),
        'init_params': result.init_params.dumps() if hasattr(result, 'init_params') else None,
        'success': result.success if hasattr(result, 'success') else None,
        'message': result.message if hasattr(result, 'message') else None,
        'nfev': result.nfev if hasattr(result, 'nfev') else None,
        'nvarys': result.nvarys if hasattr(result, 'nvarys') else None,
        'ndata': result.ndata if hasattr(result, 'ndata') else None,
        'nfree': result.nfree if hasattr(result, 'nfree') else None,
        'chisqr': result.chisqr if hasattr(result, 'chisqr') else None,
        'redchi': result.redchi if hasattr(result, 'redchi') else None,
        'aic': result.aic if hasattr(result, 'aic') else None,
        'bic': result.bic if hasattr(result, 'bic') else None,
    }

    # Save the fit result
    with open(filename, 'w') as f:
        json.dump(state, f, indent=2)

    # Save the model if requested
    if model_filename != '':
        if model_filename is None:
            model_filename = filename.replace('.json', '_model.json')
            if model_filename == filename:
                model_filename = filename.replace('.json', '') + '_model.json'

        if hasattr(result, 'model') and isinstance(result.model, TransmissionModel):
            result.model.save(model_filename)


def load_result(filename: str, model_filename: str = None, model: TransmissionModel = None):
    """
    Load a ModelResult from JSON file(s).

    This function reconstructs a fit result from saved files. It can either
    load the model from a separate file or use a provided model instance.

    Parameters
    ----------
    filename : str
        Path to the fit result JSON file.
    model_filename : str, optional
        Path to the model configuration file. If None, looks for
        filename.replace('.json', '_model.json').
    model : TransmissionModel, optional
        If provided, uses this model instead of loading from file.
        Useful when you already have the model instance.

    Returns
    -------
    dict
        A dictionary containing:
        - 'params': lmfit.Parameters with fitted values
        - 'init_params': lmfit.Parameters with initial values
        - 'model': TransmissionModel (if loaded or provided)
        - 'statistics': dict with fit statistics (chisqr, redchi, etc.)
        - All other fit result attributes

    Notes
    -----
    This function returns a dictionary instead of a full ModelResult object
    because reconstructing the complete ModelResult requires re-running the fit.
    The returned dictionary contains all the essential information from the fit.

    Examples
    --------
    >>> # Load with model
    >>> result_data = load_result('my_fit.json')
    >>> print(result_data['params'])
    >>> print(result_data['statistics']['redchi'])
    >>>
    >>> # Use the loaded model for a new fit
    >>> model = result_data['model']
    >>> new_result = model.fit(new_data, params=result_data['params'])
    """
    import json

    # Load fit result
    with open(filename, 'r') as f:
        state = json.load(f)

    # Verify version
    if state.get('version') != '1.0':
        warnings.warn(f"Loading result saved with version {state.get('version')}, "
                     f"current version is 1.0. Compatibility issues may occur.")

    # Load parameters
    params = lmfit.Parameters()
    params.loads(state['params'])

    init_params = None
    if state['init_params'] is not None:
        init_params = lmfit.Parameters()
        init_params.loads(state['init_params'])

    # Load or use provided model
    loaded_model = model
    if model is None:
        if model_filename is None:
            model_filename = filename.replace('.json', '_model.json')
            if model_filename == filename:
                model_filename = filename.replace('.json', '') + '_model.json'

        try:
            loaded_model = TransmissionModel.load(model_filename)
        except FileNotFoundError:
            warnings.warn(f"Model file {model_filename} not found. "
                         f"Returning results without model.")

    # Prepare return dictionary with all information
    result_dict = {
        'params': params,
        'init_params': init_params,
        'model': loaded_model,
        'statistics': {
            'success': state.get('success'),
            'message': state.get('message'),
            'nfev': state.get('nfev'),
            'nvarys': state.get('nvarys'),
            'ndata': state.get('ndata'),
            'nfree': state.get('nfree'),
            'chisqr': state.get('chisqr'),
            'redchi': state.get('redchi'),
            'aic': state.get('aic'),
            'bic': state.get('bic'),
        },
        'version': state.get('version'),
        'success': state.get('success'),
        'message': state.get('message'),
        'nfev': state.get('nfev'),
        'nvarys': state.get('nvarys'),
        'ndata': state.get('ndata'),
        'nfree': state.get('nfree'),
        'chisqr': state.get('chisqr'),
        'redchi': state.get('redchi'),
        'aic': state.get('aic'),
        'bic': state.get('bic'),
    }

    return result_dict