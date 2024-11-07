import numpy as np
import pandas as pd
from scipy.stats import exponnorm
import matplotlib.pyplot as plt
from scipy.special import erfc
import lmfit
import inspect

class Response:
    def __init__(self, kind="jorgensen", vary: bool = False,
                 wlstep=0.1):
        """
        Initializes the Response object with specified parameters.

        Parameters:
        kind (str): The type of response function to use. Options are 'expo_gauss' or 'none'.
        vary (bool): If True, the parameters can vary during fitting. Default is False.
        eps (float): The threshold for cutting the response array symmetrically. Default is 1.0e-6.
        tstep (float): The time step for the response function. Default is 1.56255e-9 seconds.
        nbins (int): The number of bins for the response function. Default is 300.
        """
        self.wlstep = wlstep
        self.params = lmfit.Parameters()
        self.Δλ = np.arange(-20, 20, self.wlstep)

        # Choose the response function
        if kind == "jorgensen":
            self.function = self.jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α1", value=3.67, min=0.001, max= 5, vary=vary)
            self.params.add(f"β1", value=3.06, min=0.001, max= 5, vary=vary)

        elif kind == "bem":
            self.function = self.bem_response
            self.params = lmfit.Parameters()
            self.params.add(f"α1", value=3.67, min=0.001, max= 5, vary=vary)
            self.params.add(f"β1", value=3.06, min=0.001, max= 5, vary=vary)

        elif kind == "none":
            self.function = self.empty_response
        else:
            raise NotImplementedError(f"Response kind '{kind}' is not supported. Use 'jorgensen' or 'none'.")

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
    
    def bem_response(self, α1, β1, **kwargs):
        from bem import peak_profile
        import warnings

        Δλ = np.arange(-20, 20, 0.1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            j = peak_profile.Jorgensen(alpha=[α1, 0.], beta=[β1, 0], sigma=[0, 0, 0]).calc_profile(Δλ, 1.)
        j = j/sum(j)
        return np.nan_to_num(j,0)

    def jorgensen_response(self, α1, β1, **kwargs):
        """
        Calculate the Jorgensen peak profile with a fixed self.Δλ range.

        This function implements the Jorgensen peak profile based on the formulas
        described in section 3.3.3 of Sven Vogel's thesis,
        "A Rietveld-approach for the analysis of neutron time-of-flight transmission data".

        Parameters:
        -----------
        α1 : float
            The first alpha parameter. Units: 1/angstrom.
        β1 : float
            The first beta parameter. Units: 1/angstrom.

        Returns:
        --------
        numpy.ndarray
            The calculated Jorgensen profile, normalized so that its sum is 1.
        """

        # Calculate parameters
        alpha = α1
        beta = β1
        sigma = 1.0  # Fixed sigma value

        # Calculate scale and other constants
        scale = alpha * beta / (2 * (alpha + beta))
        sigma2 = sigma * sigma
        sqrt2 = np.sqrt(2)

        # Calculate u, v, y, and z
        u = alpha / 2 * (alpha * sigma2 + 2 * self.Δλ)
        v = beta / 2 * (beta * sigma2 - 2 * self.Δλ)
        y = (alpha * sigma2 + self.Δλ) / (sqrt2 * sigma)
        z = (beta * sigma2 - self.Δλ) / (sqrt2 * sigma)

        # Calculate the profile
        term1 = np.exp(u) * erfc(y)
        term2 = np.exp(v) * erfc(z)
        
        # Avoid division by zero warnings
        term1 = np.where(np.isfinite(term1), term1, 0)
        term2 = np.where(np.isfinite(term2), term2, 0)

        profile = scale * (term1 + term2)

        # Normalize the profile
        normalized_profile = profile / np.sum(profile)

        return np.nan_to_num(normalized_profile,0.)

    def plot(self, params=None, **kwargs):
        """
        Plots the response function.

        Parameters:
        params (dict): Parameters for the response function.
        **kwargs: Additional arguments for plot customization.
        """
        ax = kwargs.pop("ax", plt.gca())
        xlabel = kwargs.pop("xlabel", "wavelength [Angstrom]")

        params = params if params else self.params
        y = self.function(**params.valuesdict())
        df = pd.Series(y, index=self.Δλ, name="Response")
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
            for i in range(3):
                self.params.add(f"b{i}", value=0., min=-1e6, max= 1e6, vary=vary)

        elif kind == "polynomial5":
            self.function = self.polynomial5_background
            for i in range(5):
                self.params.add(f"b{i}", value=0., min=-1e6, max= 1e6, vary=vary)

        elif kind == "constant":
            self.function = self.constant_background
            self.params.add('b0', value=0.0, min=-1e6,max=1e6,vary=vary)
        elif kind == "none":
            self.function = self.empty_background
        else:
            raise NotImplementedError(f"Background kind '{kind}' is not supported. Use 'none', 'constant', 'polynomial3', or 'polynomial5'.")

    def empty_background(self, wl, **kwargs):
        """
        Returns a zero background array.
        """
        return np.zeros_like(wl)

    def constant_background(self, wl, b0=0., **kwargs):
        """
        Generates a constant background.

        Parameters:
        wl (np.ndarray): Wavelength values.
        b0 (float): Constant background value.
        """
        bg = np.full_like(wl, b0)
        return np.where(bg>0,bg,0.)

    def polynomial3_background(self, wl, b0=0., b1=1., b2=0., **kwargs):
        """
        Computwls a third-degree polynomial background.

        Parameters:
        wl (np.ndarray): Energy values.
        b0 (float): Constant term.
        b1 (float): Linear term.
        b2 (float): Quadratic term.
        """
        bg = b0 + b1 * np.sqrt(wl) + b2 / np.sqrt(wl)
        return np.where(bg>0,bg,0.)

    def polynomial5_background(self, wl, b0=0., b1=1., b2=0., b3=0., b4=0., **kwargs):
        """
        Computes a fifth-degree polynomial background.

        Parameters:
        wl (np.ndarray): Wavelegth values.
        b0 (float): Constant term.
        b1 (float): Linear term.
        b2 (float): Quadratic term.
        b3 (float): Cubic term.
        b4 (float): Quartic term.
        """
        bg = b0 + b1 * np.sqrt(wl) + b2 / np.sqrt(wl) + b3 * wl + b4 * wl**2
        return np.where(bg>0,bg,0.)

    def plot(self, wl, params=None, **kwargs):
        """
        Plots the background function.

        Parameters:
        wl (np.ndarray): Wavelength values.
        params (dict): Parameters for the background function.
        """
        ax = kwargs.pop("ax", plt.gca())
        ls = kwargs.pop("ls", "--")
        color = kwargs.pop("color", "0.5")
        params = params if params else self.params
        y = self.function(wl, **params.valuesdict())
        df = pd.Series(y, index=wl, name="Background")
        df.plot(ax=ax, color=color, ls=ls, **kwargs)
