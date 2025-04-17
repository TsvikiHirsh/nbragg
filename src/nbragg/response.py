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
                 wlstep:float =0.1,tstep:float =10e-6, cut_threshold: float=1e-6):
        """
        Initializes the Response object with specified parameters.

        Parameters:
        kind (str): The type of response function to use. Options are 'jorgensen', 'square', 'square_jorgensen' or 'none'.
        vary (bool): If True, the parameters can vary during fitting. Default is False.
        wlstep (float, optional): wavelength step size (Angstrom)
        tstep (float, optional): time step size (s)
        cut_threshold (float, optional): threshold to cut the response function
        """
        self.wlstep = wlstep
        self.params = lmfit.Parameters()
        self.Δλ = np.arange(-15, 15.00001, self.wlstep)
        self.tgrid = np.arange(-0.004,0.0040001,tstep) # grid for time based response -5ms to 5ms with 10usec bin size
        self.kind = kind
        self.cut_threshold = cut_threshold

        # Choose the response function
        if kind == "jorgensen":
            self.function = self.jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max= 10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max= 10000, vary=vary)

        elif kind == "square":
            self.function = self.square_response
            self.params = lmfit.Parameters()
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max= 5000, vary=vary)

        elif kind == "square_jorgensen":
            self.function = self.square_jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max= 10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max= 10000, vary=vary)
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max= 5000, vary=vary)

        elif kind == "full_jorgensen":
            self.function = self.full_jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max= 10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max= 10000, vary=vary)
            self.params.add(f"σ1", value=2.5e-3, min=1e-6, max= 100e-3, vary=vary)

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

        response = arr[left_bound:right_bound]

        return response / np.sum(response)

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
        tof_response = self.cut_array_symmetric(tof_response,self.cut_threshold)
        return tof_response

    def jorgensen_response(self, α0, β0, σ1=None, **kwargs):
        """
        Calculates the Jorgensen peak profile function with Greek Unicode parameters.
        
        Parameters:
        -----------
        α0 : float or list/tuple
            Alpha parameter [α0₁, α0₂] or single value α0₁ (α0₂ defaults to 0)
        β0 : float or list/tuple
            Beta parameter [β0₁, β0₂] or single value β0₁ (β0₂ defaults to 0)
        σ1 : list/tuple, optional
            Sigma parameters [σ₁, σ₂, σ₃], defaults to [0, 0, 0]
        
        Returns:
        --------
        numpy.ndarray
            Normalized profile values with NaN values replaced by 0
        """
        # Handle input parameters
        if σ1 is None:
            σ1 = [0, 0, 0]
        
        # Convert single values to lists with 0 as second element
        if not isinstance(α0, (list, tuple)):
            α0 = [α0, 0]
        if not isinstance(β0, (list, tuple)):
            β0 = [β0, 0]
        
        # Generate x values
        x = self.Δλ
        
        # Calculate parameters using d-spacing of 1.0 as in original code
        d = 1.0
        α0_calc = α0[0] + α0[1]/d
        β0_calc = β0[0] + β0[1]/d**4
        σ1_calc = np.sqrt(σ1[0]**2 + (σ1[1]*d)**2 + (σ1[2]*d*d)**2)
        
        # Constants
        sqrt2 = np.sqrt(2)
        σ2 = σ1_calc * σ1_calc
        
        # Scaling factor
        scale = α0_calc * β0_calc / 2 / (α0_calc + β0_calc)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Calculate intermediate terms
            u = α0_calc/2 * (α0_calc*σ2 + 2*x)
            v = β0_calc/2 * (β0_calc*σ2 - 2*x)
            y = (α0_calc*σ2 + x)/(sqrt2*σ1_calc)
            z = (β0_calc*σ2 - x)/(sqrt2*σ1_calc)
            
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
            profile = np.nan_to_num(profile, 0)
            return profile
            # cut symmetric
            return self.cut_array_symmetric(profile,self.cut_threshold)
        
    def full_jorgensen_response(self, wl=4., α0=3.67, β0=3.04, σ1=0.6e-3, **kwargs):
        """
        Calculates the Jorgensen peak profile function.
        
        Parameters:
        -----------
        wl : array of wavelengths
        α0 : float or list/tuple
            Alpha parameter [α0, α0] or single value α0 (α0 defaults to 0)
        β0 : float or list/tuple
            Beta parameter [β0, β0] or single value β0 (β0 defaults to 0)
        σ1 : list/tuple, optional
            Sigma parameters [σ1₁, σ1₂, σ1₃], defaults to [0, 1, 0]
        
        Returns:
        --------
        numpy.ndarray
            Normalized profile values with NaN values replaced by 0
        """
        # Handle input parameters
        if not isinstance(σ1, (list, tuple)):
            σ1 = [0, σ1, 0]
        
        # Convert single values to lists with 0 as second element
        if not isinstance(α0, (list, tuple)):
            α0 = [α0, 0]
        if not isinstance(β0, (list, tuple)):
            β0 = [β0, 0]
        
        # Generate x values
        x = self.Δλ
        
        # Calculate parameters using d-spacing =wl/2
        d = wl/2.
        α0_calc = α0[0] + α0[1]/d
        β0_calc = β0[0] + β0[1]/d**4
        σ1_calc = np.sqrt(σ1[0]**2 + (σ1[1]*d)**2 + (σ1[2]*d*d)**2)
        
        # Constants
        sqrt2 = np.sqrt(2)
        σ2 = σ1_calc * σ1_calc
        
        # Scaling factor
        scale = α0_calc * β0_calc / 2 / (α0_calc + β0_calc)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Calculate intermediate terms
            u = α0_calc/2 * (α0_calc*σ2 + 2*x)
            v = β0_calc/2 * (β0_calc*σ2 - 2*x)
            y = (α0_calc*σ2 + x)/(sqrt2*σ1_calc)
            z = (β0_calc*σ2 - x)/(sqrt2*σ1_calc)
            
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
            profile = np.nan_to_num(profile, 0)

            # cut symmetric
            return self.cut_array_symmetric(profile,self.cut_threshold)
        
    def square_jorgensen_response(self, α0=3.67, β0=3, σ1=None, width=None, **kwargs):
        """
        Calculates the Jorgensen peak profile function with an additional square width parameter.
        
        Parameters:
        -----------
        α0 : float or list/tuple
            Alpha parameter [α0, α0] or single value α0 (α0 defaults to 0)
        β0 : float or list/tuple
            Beta parameter [β0, β0] or single value β0 (β0 defaults to 0)
        σ1 : list/tuple, optional
            Sigma parameters [σ₁, σ₂, σ₃], defaults to [0, 0, 0]
        width : float, optional
            Square width to broaden the response [usec]
        
        Returns:
        --------
        numpy.ndarray
            Normalized profile values with NaN values replaced by 0
        """
        # Handle input parameters
        if σ1 is None:
            σ1 = [0, 0, 0]
        
        # Convert single values to lists with 0 as second element
        if not isinstance(α0, (list, tuple)):
            α0 = [α0, 0]
        if not isinstance(β0, (list, tuple)):
            β0 = [β0, 0]
        
        # Generate x values
        x = self.Δλ
        
        # Calculate parameters using d-spacing of 1.0 as in original code
        d = 1.0
        α0_calc = α0[0] + α0[1]/d
        β0_calc = β0[0] + β0[1]/d**4
        σ1_calc = np.sqrt(σ1[0]**2 + (σ1[1]*d)**2 + (σ1[2]*d*d)**2)
        
        # Constants
        sqrt2 = np.sqrt(2)
        σ2 = σ1_calc * σ1_calc
        
        # Scaling factor
        scale = α0_calc * β0_calc / 2 / (α0_calc + β0_calc)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # Calculate intermediate terms
            u = α0_calc/2 * (α0_calc*σ2 + 2*x)
            v = β0_calc/2 * (β0_calc*σ2 - 2*x)
            y = (α0_calc*σ2 + x)/(sqrt2*σ1_calc)
            z = (β0_calc*σ2 - x)/(sqrt2*σ1_calc)
            
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
            profile = np.nan_to_num(profile, 0)

            # cut symmetric
            return self.cut_array_symmetric(profile,self.cut_threshold)
    

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
