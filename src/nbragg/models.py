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


class TransmissionModel(lmfit.Model):
    def __init__(self, cross_section, 
                        response: str = "bem",
                        background: str = "polynomial3",
                        tof_length: float = 9,
                        vary_weights: bool = None, 
                        vary_background: bool = None, 
                        vary_tof: bool = None,
                        vary_response: bool = None,
                        vary_orientation: bool = None,
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
        kwargs : dict, optional
            Additional keyword arguments for model and background parameters.

        Notes
        -----
        This model calculates the transmission function as a combination of 
        cross-section, response function, and background.
        """
        super().__init__(self.transmission, **kwargs)

        self.cross_section = cross_section
        self.tof_length = tof_length

        self.params = self._make_basic_params()
        self.params += self._make_temperature_params() # add temperature params to fit
        if vary_weights is not None:
            self.params += self._make_weight_params(vary=vary_weights)
        if vary_tof is not None:
            self.params += self._make_tof_params(vary=vary_tof,**kwargs)

        self.response = None
        if vary_response is not None:
            self.response = Response(kind=response,vary=vary_response)
            self.params += self.response.params

        self.background = None
        if vary_background is not None:
            self.background = Background(kind=background,vary=vary_background)
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
            bg = self.background.function(wl,**kwargs)
        else:
            bg = 0.

        n = self.atomic_density

        # Transmission function
        
        xs = self.cross_section(wl,**kwargs)

        if self.response != None:
            response = self.response.function(**kwargs)
            xs = convolve1d(xs,response,0)

        T = norm * np.exp(- xs * thickness * n) * (1 - bg) + bg
        return T

    def fit(self, data, params=None, wlmin=1., wlmax=6., **kwargs):
        """
        Fit the model to the data.

        Parameters
        ----------
        data : pandas.DataFrame or nbragg.data.Data
            The data to fit the model to.
        params : lmfit.Parameters, optional
            Initial parameter values for the fit. If None, the current model parameters will be used.
        wlmin : float, optional
            The minimum wavelength for fitting, by default 0.5e6.
        wlmax : float, optional
            The maximum wavelength for fitting, by default 20.e6.
        kwargs : dict, optional
            Additional arguments passed to the lmfit.Model.fit method.

        Returns
        -------
        lmfit.model.ModelResult
            The result of the fit.

        Notes
        -----
        This function applies wavelength filtering to the input data based on `wlmin` and `wlmax`,
        then fits the transmission model to the filtered data.
        """
        # self.cross_section.set_wavelength_range(wlmin,wlmax)
        if isinstance(data,pandas.DataFrame):
            data = data.query(f"{wlmin}<wavelength<{wlmax}")
            weights = kwargs.get("weights",1./data["err"].values)
            fit_result = super().fit(data["trans"].values, params=params or self.params, weights=weights, wl=data["wavelength"].values, **kwargs)
        elif isinstance(data,Data):
            data = data.table.query(f"{wlmin}<wavelength<{wlmax}")
            weights = kwargs.get("weights",1./data["err"].values)
            fit_result = super().fit(data["trans"].values, params=params or self.params, weights=weights, wl=data["wavelength"].values, **kwargs)
        else:
            # Perform the fit using the parent class's fit method
            fit_result = super().fit(data, params=params or self.params, **kwargs)
        self.fit_result = fit_result
        # switch method names
        # fit_result.plot_results = deepcopy(fit_result.plot)
        fit_result.plot = self.plot

        # return TransmissionModelResult(fit_result, params or self.params)
        return fit_result
    
    def plot(self, plot_bg: bool = True, **kwargs):
        """
        Plot the results of the fit.

        Parameters
        ----------
        plot_bg : bool, optional
            Whether to include the background in the plot, by default True.
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
        """
        fig, ax = plt.subplots(2,1,sharex=True,height_ratios=[3.5,1],figsize=(6,5))
        wavelength = self.fit_result.userkws["wl"]
        data = self.fit_result.data
        err = 1./self.fit_result.weights
        best_fit = self.fit_result.best_fit
        residual = self.fit_result.residual
        color = kwargs.pop("color","seagreen")
        ecolor = kwargs.pop("ecolor","0.8")
        ms = kwargs.pop("ms",2)
        ax[0].errorbar(wavelength,data,err,marker="o",color=color,ms=ms,zorder=-1,ecolor=ecolor,label="Best fit")  
        ax[0].plot(wavelength,best_fit,color="0.2",label="Data") 
        ax[0].set_ylabel("Transmission")
        ax[0].set_title(self.cross_section.name)
        ax[1].plot(wavelength,residual,color=color)
        ax[1].set_ylabel("Residuals [1σ]")
        ax[1].set_xlabel("λ [Å]")
        if plot_bg and self.background.params:
            self.background.plot(wl=wavelength,ax=ax[0],params=self.fit_result.params,**kwargs)
            ax[0].legend(["Best fit","Background","Data"], fontsize=9,reverse=True,title=f"χ$^2$: {self.fit_result.redchi:.2f}")
        else:
            ax[0].legend(["Best fit","Data"], fontsize=9,reverse=True,title=f"χ$^2$: {self.fit_result.redchi:.2f}")
        plt.subplots_adjust(hspace=0.05)
        
        return ax
    
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
        for i,material in enumerate(self.cross_section.materials):
            mos = self.cross_section.materials[material].get("mos",None)
            if mos: 
                params.add(f"η{i+1}", value=mos, min=0.001, max= 50, vary=vary)
                theta = self.cross_section.materials[material].get("theta",0.)
                params.add(f"θ{i+1}", value=theta, min=0., max= 180, vary=vary)
                phi = self.cross_section.materials[material].get("phi",0.)
                params.add(f"ϕ{i+1}", value=phi, min=0, max= 360, vary=vary)
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
        for i,material in enumerate(self.cross_section.materials):
            temp = self.cross_section.materials[material].get("temp",None)
            if temp: 
                params.add(f"temp", value=temp, min=77., max= 1000, vary=vary)
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
        params.add(f"thickness", value=1., min=0., max= 5., vary=vary)
        params.add(f"norm", value=1., min=0.1, max= 10., vary=vary)
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
        params.add("L0", value=L0, min=0.5, max= 1.5, vary=vary)
        params.add("t0", value=t0, min=-5e-6, max=5e6,vary=vary)
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
        # weight_series = deepcopy(self.cross_section.weights)
        # weight_series.index = weight_series.index.str.replace("-", "")
        # param_names = weight_series.index
        # N = len(weight_series)

        # Normalize the input weights to sum to 1
        # weights = np.array(weight_series / weight_series.sum(), dtype=np.float64)


        weights = np.array([self.cross_section.materials[phase]["weight"] for phase in self.cross_section.materials])
        param_names = [phase.replace("-", "") for phase in self.cross_section.materials]

        N = len(weights)
        if N == 1:
            # Special case: if N=1, the weight is always 1
            params.add(f'{param_names[0]}', value=1., vary=False)
        else:

            last_weight = weights[-1]
            # Add (N-1) free parameters corresponding to the first (N-1) items
            for i, name in enumerate(param_names[:-1]):
                initial_value = weights[i]  # Use weight values
                params.add(f'p{i+1}',value=np.log(weights[i]/last_weight),min=-14,max=14) # limit to 1ppm
            
            # Define the normalization expression
            normalization_expr = ' + '.join([f'exp(p{i+1})' for i in range(N-1)])
            
            # Add weights based on the free parameters
            for i, name in enumerate(param_names[:-1]):
                params.add(f'{name}', expr=f'exp(p{i+1}) / (1 + {normalization_expr})')
            
            # The last weight is 1 minus the sum of the previous weights
            params.add(f'{param_names[-1]}', expr=f'1 / (1 + {normalization_expr})')

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

