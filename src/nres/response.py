import numpy as np
import pandas as pd
from scipy.stats import exponnorm
import matplotlib.pyplot as plt
import lmfit
import inspect

class Response:
    def __init__(self, kind="expo_gauss", vary: bool = False, eps: float = 1.0e-6,
                 tstep=1.56255e-9, nbins=300):
        """
        Initializes the Response object with specified parameters.

        Parameters:
        kind (str): The type of response function to use. Options are 'expo_gauss' or 'none'.
        vary (bool): If True, the parameters can vary during fitting. Default is False.
        eps (float): The threshold for cutting the response array symmetrically. Default is 1.0e-6.
        tstep (float): The time step for the response function. Default is 1.56255e-9 seconds.
        nbins (int): The number of bins for the response function. Default is 300.
        """
        self.tstep = tstep
        self.grid = np.arange(-nbins, nbins + 1, 1)
        self.tgrid = self.grid * self.tstep
        self.eps = eps
        self.params = lmfit.Parameters()

        # Choose the response function
        if kind == "expo_gauss":
            self.function = self.expogauss_response
            self.params = lmfit.Parameters()
            self.params.add_many(
                ('K', 1.0, True, 0.0001),  # Exponential shape parameter
                ('x0', 1e-9, vary),         # Location parameter (Gaussian)
                ('τ', 1e-9, True, 1e-10)    # Exponential scale parameter
            )
        elif kind == "none":
            self.function = self.empty_response
        else:
            raise NotImplementedError(f"Response kind '{kind}' is not supported. Use 'expo_gauss' or 'none'.")

    def register_response(self, response_func, lmfit_params=None, **kwargs):
        """
        Registers a new response using any scipy.stats function.

        Parameters:
        response_func (function): A function from scipy.stats, e.g., exponnorm.pdf.
        lmfit_params (lmfit.Parameters): Optional lmfit.Parameters to define limits and vary.
        kwargs: Default parameter values for the response function.
        """
        

        # Detect parameters of the response function
        sig_params = inspect.signature(response_func).parameters
        for param, default in kwargs.items():
            if param in sig_params:
                self.params.add(param, value=default, vary=True)
            else:
                raise ValueError(f"Parameter '{param}' not found in the response function signature.")
            
        self.function = response_func.pdf(self.tgrid)

        # Use optional lmfit.Parameters to customize limits and vary
        if lmfit_params:
            for name, param in lmfit_params.items():
                if name in self.params:
                    self.params[name].set(value=param.value, vary=param.vary, min=param.min, max=param.max)

    def cut_array_symmetric(self, arr, threshold):
        """
        Symmetrically cuts the array based on a threshold.

        Parameters:
        arr (np.ndarray): Input array to be cut.
        threshold (float): The threshold value for cutting the array.

        Returns:
        np.ndarray: Symmetrically cut array with an odd number of elements.
        """
        if len(arr) % 2 == 0:
            raise ValueError("Input array length must be odd.")

        center_idx = len(arr) // 2
        left_idx = np.argmax(arr[:center_idx][::-1] < threshold)
        right_idx = np.argmax(arr[center_idx:] < threshold)
        
        left_bound = center_idx - max(left_idx, right_idx)
        right_bound = center_idx + max(left_idx, right_idx) + 1  # Ensure odd length

        return arr[left_bound:right_bound]

    def empty_response(self, **kwargs):
        """
        Returns an empty response [0.0, 1.0, 0.0].
        """
        return np.array([0., 1., 0.])

    def expogauss_response(self, K=0.01, x0=0., τ=1.0e-9, **kwargs):
        """
        Computes the exponential-Gaussian response function.

        Parameters:
        K (float): Shape parameter for the exponential.
        x0 (float): Location parameter for the Gaussian.
        τ (float): Scale parameter for the exponential.

        Returns:
        np.ndarray: Normalized response array.
        """
        response = exponnorm.pdf(self.tgrid, K, loc=x0, scale=τ)
        response /= np.sum(response)
        return self.cut_array_symmetric(response, self.eps)

    def plot(self, params=None, **kwargs):
        """
        Plots the response function.

        Parameters:
        params (dict): Parameters for the response function.
        **kwargs: Additional arguments for plot customization.
        """
        ax = kwargs.pop("ax", plt.gca())
        xlabel = kwargs.pop("xlabel", "t [sec]")

        params = params if params else self.params
        y = self.function(**params.valuesdict())
        tof = np.arange(-len(y) // 2 + 1, len(y) // 2 + 1) * self.tstep
        df = pd.Series(y, index=tof, name="Response")
        df.plot(ax=ax, xlabel=xlabel, **kwargs)


class Background:
    def __init__(self, kind="expo_norm", vary: bool = False):
        """
        Initializes the Background object with specified parameters.

        Parameters:
        kind (str): Type of background function ('constant', 'polynomial3', 'polynomial5', or 'none').
        vary (bool): If True, the parameters can vary during fitting.
        """
        self.params = lmfit.Parameters()
        if kind == "polynomial3":
            self.function = self.polynomial3_background
            self.params.add_many(
                ('b0', 0.0, vary),
                ('b1', 0.0, vary),
                ('b2', 0.0, vary)
            )
        elif kind == "polynomial5":
            self.function = self.polynomial5_background
            self.params.add_many(
                ('b0', 0.0, vary),
                ('b1', 0.0, vary),
                ('b2', 0.0, vary),
                ('b3', 0.0, vary),
                ('b4', 0.0, vary)
            )
        elif kind == "constant":
            self.function = self.constant_background
            self.params.add('b0', 0.0, vary=vary)
        elif kind == "none":
            self.function = self.empty_background
        else:
            raise NotImplementedError(f"Background kind '{kind}' is not supported. Use 'none', 'constant', 'polynomial3', or 'polynomial5'.")

    def empty_background(self, E, **kwargs):
        """
        Returns a zero background array.
        """
        return np.zeros_like(E)

    def constant_background(self, E, b0=0., **kwargs):
        """
        Generates a constant background.

        Parameters:
        E (np.ndarray): Energy values.
        b0 (float): Constant background value.
        """
        return np.full_like(E, b0)

    def polynomial3_background(self, E, b0=0., b1=1., b2=0., **kwargs):
        """
        Computes a third-degree polynomial background.

        Parameters:
        E (np.ndarray): Energy values.
        b0 (float): Constant term.
        b1 (float): Linear term.
        b2 (float): Quadratic term.
        """
        return b0 + b1 * np.sqrt(E) + b2 / np.sqrt(E)

    def polynomial5_background(self, E, b0=0., b1=1., b2=0., b3=0., b4=0., **kwargs):
        """
        Computes a fifth-degree polynomial background.

        Parameters:
        E (np.ndarray): Energy values.
        b0 (float): Constant term.
        b1 (float): Linear term.
        b2 (float): Quadratic term.
        b3 (float): Cubic term.
        b4 (float): Quartic term.
        """
        return b0 + b1 * np.sqrt(E) + b2 / np.sqrt(E) + b3 * E + b4 * E**2

    def plot(self, E, params=None, **kwargs):
        """
        Plots the background function.

        Parameters:
        E (np.ndarray): Energy values.
        params (dict): Parameters for the background function.
        """
        ax = kwargs.pop("ax", plt.gca())
        ls = kwargs.pop("ls", "--")
        color = kwargs.pop("color", "0.5")
        params = params if params else self.params
        y = self.function(E, **params.valuesdict())
        df = pd.Series(y, index=E, name="Background")
        df.plot(ax=ax, color=color, ls=ls, **kwargs)
