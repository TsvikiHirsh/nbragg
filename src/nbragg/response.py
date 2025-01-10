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
from functools import wraps

class Response:
    def __init__(self, kind="jorgensen", vary: bool = False,
                 tstep: float = 10e-6, λstep:float = None, flight_path_length: float = 9.,
                 cut_symmetric: bool = False, cut_threshold: float = 1e-6):
        """
        Initializes the Response object with specified parameters.

        Parameters:
        kind (str): The type of response function to use. Options are 'jorgensen', 'square', 'square_jorgensen' or 'none'.
        vary (bool): If True, the parameters can vary during fitting. Default is False.
        cut_symmetric (bool): If True, applies symmetric cutting to the response function output.
        cut_threshold (float): Threshold value for symmetric cutting.
        """
        self.params = lmfit.Parameters()
        if λstep != None:
            self.λgrid = np.arange(-2, 2, λstep)
            self.tgrid = self.λgrid / 3956.034*flight_path_length
        else:
            self.tgrid = np.arange(-0.005, 0.005, tstep)
            self.λgrid = self.tgrid * 3956.034/flight_path_length
        self.kind = kind
        self.cut_symmetric = cut_symmetric
        self.cut_threshold = cut_threshold

        def symmetric_cut_decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                result = func(*args, **kwargs)
                if self.cut_symmetric:
                    return self.cut_array_symmetric(result, self.cut_threshold)
                return result
            return wrapper

        # Choose the response function and apply decorator if needed
        if kind == "jorgensen":
            self.function = symmetric_cut_decorator(self.jorgensen_response)
            self.params = lmfit.Parameters()
            self.params.add(f"α1", value=3.67, min=0.001, max=1000, vary=vary)
            self.params.add(f"β1", value=3.06, min=0.001, max=1000, vary=vary)

        elif kind == "square":
            self.function = symmetric_cut_decorator(self.square_response)
            self.params = lmfit.Parameters()
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max=5000, vary=vary)

        elif kind == "square_jorgensen":
            self.function = symmetric_cut_decorator(self.square_jorgensen_response)
            self.params = lmfit.Parameters()
            self.params.add(f"α1", value=3.67, min=0.001, max=1000, vary=vary)
            self.params.add(f"β1", value=3.06, min=0.001, max=1000, vary=vary)
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max=5000, vary=vary)

        elif kind == "none":
            self.function = symmetric_cut_decorator(self.empty_response)
        else:
            raise NotImplementedError(f"Response kind '{kind}' is not supported. Use 'jorgensen' or 'square_jorgensen', or 'square' or 'none'.")

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
            
        def symmetric_cut_decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                result = func(*args, **kwargs)
                if self.cut_symmetric:
                    return self.cut_array_symmetric(result, self.cut_threshold)
                return result
            return wrapper

        self.function = symmetric_cut_decorator(lambda *args, **kwargs: response_func.pdf(self.tgrid))

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
        # Find center index
        center_idx = len(arr) // 2
        
        # Find indices where values fall below threshold
        below_threshold = arr < threshold
        
        # Find the first index from left and right where value falls below threshold
        left_indices = np.where(below_threshold[:center_idx])[0]
        right_indices = np.where(below_threshold[center_idx:])[0] + center_idx
        
        if len(left_indices) == 0 or len(right_indices) == 0:
            return arr
        
        # Get the leftmost and rightmost indices
        left_bound = left_indices[-1]
        right_bound = right_indices[0]
        
        # Ensure symmetry by taking the smaller distance to center
        distance_to_center = min(center_idx - left_bound, right_bound - center_idx)
        
        # Calculate new bounds
        new_left = center_idx - distance_to_center
        new_right = center_idx + distance_to_center + 1  # +1 to include the right bound
        
        return arr[new_left:new_right]


    def empty_response(self, **kwargs):
        """
        Returns an empty response [0.0, 1.0, 0.0].
        """
        return np.array([0., 1., 0.])

    def square_response(self, width=10, center: bool = True, **kwargs):
        """
        Returns a square response in time with a given width [usec]
        """
        width = width*1e-6 # convert to sec
        loc = -0.5*width if center else 0.
        tof_response = uniform.pdf(self.tgrid,loc=loc,scale=width)
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
        x = self.λgrid
        
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
        
    def square_jorgensen_response(self, α1=3.67, β1=3, σ=None, width=None, **kwargs):
        """
        Calculates the Jorgensen peak profile function with an additional square width parameter.
        
        Parameters:
        -----------
        α1 : float or list/tuple
            Alpha parameter [α1₁, α1₂] or single value α1₁ (α1₂ defaults to 0)
        β1 : float or list/tuple
            Beta parameter [β1₁, β1₂] or single value β1₁ (β1₂ defaults to 0)
        σ : list/tuple, optional
            Sigma parameters [σ₁, σ₂, σ₃], defaults to [0, 0, 0]
        width : float, optional
            Square width to broaden the response [usec]
        
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
        x = self.λgrid
        
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
            
            # Apply square width broadening if specified
            if width is not None:
                # Convert width to seconds
                width_sec = width * 1e-6
                
                # Create square broadening function
                loc = -0.5 * width_sec
                square_broadening = uniform.pdf(self.tgrid, loc=loc, scale=width_sec)
                
                # Convolve Jorgensen profile with square broadening
                from scipy.ndimage import convolve1d
                profile = convolve1d(profile, square_broadening, 0)
            
            # Normalize
            profile = profile / np.sum(profile)
            
            # Replace any NaN values with 0
            return np.nan_to_num(profile, 0)
    

    def plot(self, params=None, by_tof=False, **kwargs):
        """
        Plots the response function.

        Parameters:
        params (dict): Parameters for the response function.
        by_tof (bool): If True plots the response function by time of flight, otherwise use wavelength bininng
        **kwargs: Additional arguments for plot customization.
        """
        ax = kwargs.pop("ax", plt.gca())
        

        params = params if params else self.params
        y = self.function(**params.valuesdict())
        if by_tof:
            xlabel = kwargs.pop("xlabel", "TOF [usec]")
            df = pd.Series(y, index=self.tgrid*1e6, name="Response")
        else:
            xlabel = kwargs.pop("xlabel", "wavelength [Angstrom]")
            df = pd.Series(y, index=self.λgrid, name="Response")
        
        df.plot(ax=ax, xlabel=xlabel, **kwargs)



class Background:
    def __init__(self, kind="expo_norm", vary: bool = False):
        """
        Initializes the Background object with specified parameters.

        Parameters:
        kind (str): Type of background function ('constant', 'polynomial3', 'polynomial5','sample_dependent' or 'none').
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
        elif kind == "sample_dependent":
            self.function = self.polynomial3_background
            self.params.add(f"k", value=1., min=1., max= 10., vary=vary)
            for i in range(3):
                self.params.add(f"b{i}", value=0., min=-1e6, max= 1e6, vary=vary)


        elif kind == "constant":
            self.function = self.constant_background
            self.params.add('b0', value=0.0, min=-1e6,max=1e6,vary=vary)
        elif kind == "none":
            self.function = self.empty_background
        else:
            raise NotImplementedError(f"Background kind '{kind}' is not supported. Use 'none', 'constant', 'polynomial3', 'sample_dependent', or 'polynomial5'.")

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
        k = params.get("k",1.) # sample dependent parameter if present
        y = k*self.function(wl, **params.valuesdict())
        df = pd.Series(y, index=wl, name="Background")
        df.plot(ax=ax, color=color, ls=ls, **kwargs)
