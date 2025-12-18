"""
Save/load functionality for TransmissionModel and fit results.

This module provides an IOMixin class that can be inherited by TransmissionModel
to add save/load capabilities, as well as module-level functions for saving and
loading fit results.
"""

import json
import os
import warnings
from typing import TYPE_CHECKING

import lmfit

from nbragg.cross_section import CrossSection
from nbragg.response import Response, Background

if TYPE_CHECKING:
    from nbragg.models import TransmissionModel


class IOMixin:
    """
    Mixin class providing save/load functionality for TransmissionModel.

    This mixin adds methods to save and load model configurations to/from JSON files,
    avoiding pickle-related issues with NCrystal objects.
    """

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


# Module-level helper functions

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
    from nbragg.models import TransmissionModel

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


def load_result(filename: str, model_filename: str = None, model: 'TransmissionModel' = None):
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
    from nbragg.models import TransmissionModel

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
