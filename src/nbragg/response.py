import numpy as np
import pandas as pd
from scipy.stats import exponnorm
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.stats import uniform
import lmfit
import inspect
import warnings
import nbragg.utils as utils

class Response:
    def __init__(self, kind="jorgensen", vary: bool = False,
                 wlstep:float =0.1,tstep:float =10e-6):
        """
        Initializes the Response object with specified parameters.

        Parameters:
        kind (str): The type of response function to use. Options are 'expo_gauss' or 'none'.
        vary (bool): If True, the parameters can vary during fitting. Default is False.
        """
        self.wlstep = wlstep
        self.params = lmfit.Parameters()
        self.Δλ = np.arange(-20, 20, self.wlstep)
        self.tgrid = np.arange(-0.005,0.005,tstep) # grid for time based response -5ms to 5ms with 10usec bin size
        self.kind = kind

        # Choose the response function
        if kind == "jorgensen":
            self.function = self.jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α1", value=3.67, min=0.001, max= 1000, vary=vary)
            self.params.add(f"β1", value=3.06, min=0.001, max= 1000, vary=vary)

        elif kind == "square":
            self.function = self.square_response
            self.params = lmfit.Parameters()
            self.params.add(f"width", value=1, min=tstep*1e6, max= 5000, vary=vary)

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

    def square_response(self, width=10, **kwargs):
        """
        Returns a square response in time with a given width [usec]
        """
        width = width*1e-6 # convert to sec
        tof_response = uniform.pdf(self.tgrid,scale=width)
        tof_response /= np.sum(tof_response)
        return tof_response

    def jorgensen_response(self, α1, β1, σ=None, **kwargs):
        """
        Calculates the Jorgensen peak profile function with Greek Unicode parameters.
        
        Parameters:
        -----------
        α1 : float or list/tuple
            Alpha parameter [α1₁, α1₂] or single value α1₁ (α1₂ defaults to 0)
        β1 : float or list/tuple
            Beta parameter [β1₁, β1₂] or single value β1₁ (β1₂ defaults to 0)
        σ : list/tuple, optional
            Sigma parameters [σ₁, σ₂, σ₃], defaults to [0, 0, 0]
        
        Returns:
        --------
        numpy.ndarray
            Normalized profile values with NaN values replaced by 0
        """
        # Handle input parameters
        if σ is None:
            σ = [0, 0, 0]
        
        # Convert single values to lists with 0 as second element
        if not isinstance(α1, (list, tuple)):
            α1 = [α1, 0]
        if not isinstance(β1, (list, tuple)):
            β1 = [β1, 0]
        
        # Generate x values
        x = self.Δλ
        
        # Calculate parameters using d-spacing of 1.0 as in original code
        d = 1.0
        α1_calc = α1[0] + α1[1]/d
        β1_calc = β1[0] + β1[1]/d**4
        σ_calc = np.sqrt(σ[0]**2 + (σ[1]*d)**2 + (σ[2]*d*d)**2)
        
        # Constants
        sqrt2 = np.sqrt(2)
        σ2 = σ_calc * σ_calc
        
        # Scaling factor
        scale = α1_calc * β1_calc / 2 / (α1_calc + β1_calc)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Calculate intermediate terms
            u = α1_calc/2 * (α1_calc*σ2 + 2*x)
            v = β1_calc/2 * (β1_calc*σ2 - 2*x)
            y = (α1_calc*σ2 + x)/(sqrt2*σ_calc)
            z = (β1_calc*σ2 - x)/(sqrt2*σ_calc)
            
            # Calculate profile with special handling for numerical stability
            term1 = np.exp(u) * erfc(y)
            term2 = np.exp(v) * erfc(z)
            
            # Zero out terms where erfc is zero to avoid NaN
            term1[erfc(y) == 0] = 0
            term2[erfc(z) == 0] = 0
            
            # Calculate profile
            profile = scale * (term1 + term2)
            
            # Normalize
            profile = profile / np.sum(profile)
            
            # Replace any NaN values with 0
            return np.nan_to_num(profile, 0)

    def plot(self, params=None, **kwargs):
        """
        Plots the response function.

        Parameters:
        params (dict): Parameters for the response function.
        **kwargs: Additional arguments for plot customization.
        """
        ax = kwargs.pop("ax", plt.gca())
        

        params = params if params else self.params
        y = self.function(**params.valuesdict())
        if self.kind == "jorgensen":
            xlabel = kwargs.pop("xlabel", "wavelength [Angstrom]")
            df = pd.Series(y, index=self.Δλ, name="Response")
            df.plot(ax=ax, xlabel=xlabel, **kwargs)
        elif self.kind == "square":
            xlabel = kwargs.pop("xlabel", "tof [usec]")
            df = pd.Series(y, index=self.tgrid, name="Response")
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
