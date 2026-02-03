"""
Grouped fit functionality for spatially-resolved neutron transmission data.

This module contains the GroupedFitResult class for storing and analyzing
fit results from grouped/spatially-resolved data, along with helper functions
for parallel fitting.
"""

import lmfit
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, TYPE_CHECKING
import json
import ast

if TYPE_CHECKING:
    from nbragg.models import TransmissionModel


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

    # Stages summary if present (for rietveld/staged fits)
    if hasattr(fit_result, 'stages_summary') and fit_result.stages_summary is not None:
        import pandas as pd
        stages_summary = fit_result.stages_summary
        # If it's a Styler, extract the underlying DataFrame
        if isinstance(stages_summary, pd.io.formats.style.Styler):
            stages_summary = stages_summary.data
        if isinstance(stages_summary, pd.DataFrame):
            result_dict['stages_summary_json'] = stages_summary.to_json(orient='split')

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

    # Restore stages summary if present
    if 'stages_summary_json' in result_dict:
        import pandas as pd
        from io import StringIO
        result.stages_summary = pd.read_json(StringIO(result_dict['stages_summary_json']), orient='split')

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
    # Import here to avoid circular dependency
    from nbragg.models import TransmissionModel

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

    def plot_stage_progression(self, index, param_name, ax=None, **kwargs):
        """
        Plot the progression of a parameter across fitting stages for a specific group.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to plot.
        param_name : str
            The name of the parameter to plot (e.g., 'norm', 'thickness', 'a').
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None, a new figure is created.
        **kwargs
            Additional keyword arguments for plotting.

        Returns:
        --------
        matplotlib.axes.Axes
            The axes containing the plot.
        """
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        fit_result = self.results[normalized_index]

        if not hasattr(fit_result, 'fit_stages') or not fit_result.fit_stages:
            raise ValueError(f"No stage results available for index {index}. "
                           "Run a multi-stage fit (method='rietveld' or 'staged') first.")

        values = []
        stderrs = []
        stage_numbers = list(range(1, len(fit_result.fit_stages) + 1))

        for stage_result in fit_result.fit_stages:
            if hasattr(stage_result, 'params') and param_name in stage_result.params:
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
        ax.set_title(f"Progression of {param_name} Across Stages (index={index})")
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        return ax

    def plot_chi2_progression(self, index, ax=None, **kwargs):
        """
        Plot the progression of reduced chi-squared across fitting stages for a specific group.

        Parameters:
        -----------
        index : int, tuple, or str
            The group index to plot.
        ax : matplotlib.axes.Axes, optional
            The axes to plot on. If None, a new figure is created.
        **kwargs
            Additional keyword arguments for plotting.

        Returns:
        --------
        matplotlib.axes.Axes
            The axes containing the plot.
        """
        normalized_index = self._normalize_index(index)

        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        fit_result = self.results[normalized_index]

        if not hasattr(fit_result, 'fit_stages') or not fit_result.fit_stages:
            raise ValueError(f"No stage results available for index {index}. "
                           "Run a multi-stage fit (method='rietveld' or 'staged') first.")

        chi2_values = []
        stage_numbers = list(range(1, len(fit_result.fit_stages) + 1))

        for stage_result in fit_result.fit_stages:
            chi2_values.append(stage_result.redchi if hasattr(stage_result, 'redchi') else np.nan)

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))

        color = kwargs.pop("color", "seagreen")
        ax.plot(stage_numbers, chi2_values, "o-", color=color, **kwargs)
        ax.set_xlabel("Stage Number")
        ax.set_ylabel("Reduced chi2")
        ax.set_title(f"Reduced chi2 Progression Across Stages (index={index})")
        ax.grid(True, linestyle="--", alpha=0.7)
        plt.tight_layout()
        return ax

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
        import pandas as pd

        normalized_index = self._normalize_index(index)
        if normalized_index not in self.results:
            raise ValueError(f"Index {index} not found. Available indices: {self.indices}")

        fit_result = self.results[normalized_index]

        # Check if it's a proper ModelResult or a SimpleNamespace
        if hasattr(fit_result, '_repr_html_'):
            html = fit_result._repr_html_()
        else:
            # If it's a SimpleNamespace (from loaded file), create a basic HTML report
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
        # Import here to avoid circular dependency
        from nbragg.models import TransmissionModel

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
        # Import here to avoid circular dependency
        from nbragg.models import TransmissionModel

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
                import pandas as pd
                from io import StringIO

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
                    result.stages_summary = pd.read_json(StringIO(result_data['stages_summary']), orient='split')

            grouped_result.add_result(idx, result)

        return grouped_result
