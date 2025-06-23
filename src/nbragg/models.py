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
import warnings
import ipywidgets as widgets
from IPython.display import display
from matplotlib.patches import Rectangle


class TransmissionModel(lmfit.Model):
    def __init__(self, cross_section, 
                        params: "lmfit.Parameters" = None,
                        response: str = "jorgensen",
                        background: str = "polynomial3",
                        tof_length: float = 9,
                        vary_weights: bool = None, 
                        vary_background: bool = None, 
                        vary_tof: bool = None,
                        vary_response: bool = None,
                        vary_orientation: bool = None,
                        vary_lattice: bool = None,
                        vary_extinction: bool = None,
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
        vary_lattice: bool, optional
            It True, allows the lattice parameters of the material to be varied 
        vary_extinction: bool, optional
            It True, allows the extinction parameters of the material to be varied (requires the CrysExtn plugin to be installed)
        kwargs : dict, optional
            Additional keyword arguments for model and background parameters.

        Notes
        -----
        This model calculates the transmission function as a combination of 
        cross-section, response function, and background.
        """
        super().__init__(self.transmission, **kwargs)

        # make a new instance of the cross section
        self.cross_section = CrossSection(cross_section,
                                          name=cross_section.name,
                                          total_weight=cross_section.total_weight)
        # update atomic density
        self.cross_section.atomic_density = cross_section.atomic_density                                          
        self._materials = self.cross_section.materials
        self.tof_length = tof_length

        if params!=None:
            self.params = params.copy()
        else:
            self.params = lmfit.Parameters()
        if "thickness" not in self.params and "norm" not in self.params:
            self.params += self._make_basic_params()
        if "temp" not in self.params:
            self.params += self._make_temperature_params() # add temperature params to fit
        if vary_weights is not None:
            self.params += self._make_weight_params(vary=vary_weights)
        if vary_tof is not None:
            self.params += self._make_tof_params(vary=vary_tof,**kwargs)
        if vary_lattice is not None:
            self.params += self._make_lattice_params(vary=vary_lattice)
        if vary_extinction is not None:
            self.params += self._make_extinction_params(vary=vary_extinction)


        self.response = None
        if vary_response is not None:
            self.response = Response(kind=response,vary=vary_response)
            if list(self.response.params.keys())[0] in self.params:
                for param_name in self.params.keys():
                    self.params[param_name].vary = vary_response 
            else:
                self.params += self.response.params


        self.background = None
        if vary_background is not None:
            self.background = Background(kind=background,vary=vary_background)
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

    def fit(self, data, params=None, wlmin:float = 1., wlmax:float = 6., 
                method:str ="leastsq",
                xtol: float = None, ftol: float = None, gtol: float = None,
                **kwargs):
        """
        Fit the model to the data.

        Parameters
        ----------
        data : pandas.DataFrame or nbragg.data.Data
            The data to fit the model to.
        params : lmfit.Parameters, optional
            Initial parameter values for the fit. If None, the current model parameters will be used.
        wlmin : float, optional
            The minimum wavelength for fitting, by default 1.0.
        wlmax : float, optional
            The maximum wavelength for fitting, by default 6.0.
        method : str, optional
            the method name, default is leastsq. For lattice variation the recommended fit method is "nelder", other valid options are "brute","cobyla","powell" or others specified in lmfit.minimize docstring
        xtol : float, optional
            Relative tolerance for changes in the parameters. The optimizer stops when the relative
            changes in parameter values are smaller than `xtol`. Default is None.
        ftol : float, optional
            Relative tolerance for the cost function (e.g., sum of squared residuals). The optimizer
            stops when the relative change in the cost function is smaller than `ftol`. Default is None.
        gtol : float, optional
            Tolerance for the gradient of the cost function with respect to the parameters. The optimizer
            stops when the gradient norm falls below `gtol`. Default is None.
        kwargs : dict, optional
            Additional arguments passed to the `lmfit.Model.fit` method.

        Returns
        -------
        lmfit.model.ModelResult
            The result of the fit, containing optimized parameters and fitting statistics.

        Notes
        -----
        - This function applies wavelength filtering to the input data based on `wlmin` and `wlmax`,
        then fits the transmission model to the filtered data.
        - The `xtol`, `ftol`, and `gtol` parameters are passed to the underlying optimization method 
        in `lmfit.Model.fit` via the `fit_kws` argument.
        If not specified, the default tolerances for the optimizer are used.
        """
        # Update fit_kws to include xtol, ftol, gtol
        fit_kws = kwargs.pop("fit_kws", {})  # Extract existing fit_kws or initialize an empty dict
        fit_kws.setdefault("xtol", xtol) if xtol is not None else None
        fit_kws.setdefault("ftol", ftol) if ftol is not None else None
        fit_kws.setdefault("gtol", gtol) if gtol is not None else None

        # Pass fit_kws back into kwargs
        kwargs["fit_kws"] = fit_kws

        # Apply wavelength filtering and weights
        if isinstance(data, pandas.DataFrame):
            data = data.query(f"{wlmin} < wavelength < {wlmax}")
            weights = kwargs.get("weights", 1. / data["err"].values)
            fit_result = super().fit(
                data["trans"].values,
                params=params or self.params,
                weights=weights,
                wl=data["wavelength"].values,
                method = method,
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
                method = method,
                **kwargs
            )

        else:
            # Perform the fit using the parent class's fit method
            fit_result = super().fit(
                data,
                params=params or self.params,
                method = method,
                **kwargs
            )

        # Store and modify fit results
        self.fit_result = fit_result
        fit_result.plot = self.plot  # Modify the plot method
        fit_result.plot_total_xs = self.plot_total_xs
        if self.response != None:
            fit_result.response = self.response
            fit_result.response.params = fit_result.params
        if self.background != None:
            fit_result.background = self.background

        return fit_result
    
    def plot(self, data=None, plot_bg: bool = True,
            plot_dspace: bool = False, dspace_min: float = 1,
            dspace_label_pos: float = 0.99, **kwargs):
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
        fig, ax = plt.subplots(2, 1, sharex=True, height_ratios=[3.5, 1], figsize=(6, 5))
        
        # Check if we have fit results or should use model
        if hasattr(self, "fit_result") and self.fit_result is not None:
            # Use fit results
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
                # You might want to define a default wavelength range here
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

    def plot_total_xs(self, plot_bg: bool = True, 
                    plot_dspace: bool = False, 
                    dspace_min: float = 1, 
                    dspace_label_pos: float = 0.99, 
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
        kwargs : dict, optional
            Additional plot settings like color, marker size, etc.

        Returns
        -------
        matplotlib.axes.Axes
            The axes of the plot.

        Notes
        -----
        This function generates a plot showing the total cross-section data, 
        the best-fit curve, and residuals. If `plot_bg` is True, it will also 
        plot the background function.
        """
        # Prepare data for cross-section calculation
        wavelength = self.fit_result.userkws["wl"]
        if "k" in self.fit_result.params:
            k = self.fit_result.params['k'].value
        else:
            k = 1.
        if plot_bg and self.background:
            bg = self.background.function(wavelength,**self.fit_result.params)
        else:
            bg = 0.
        norm = self.fit_result.params['norm'].value
        n = self.atomic_density
        thickness = self.fit_result.params['thickness'].value

        # Calculate cross-section data
        data_xs = -1. / n / thickness * np.log((self.fit_result.data - k * bg) / norm / (1. - bg))

        fig, ax = plt.subplots(2, 1, sharex=True, height_ratios=[3.5, 1], figsize=(6, 5))

        
        # Calculate best fit and residuals for cross-section
        xs = self.cross_section(wavelength,**self.fit_result.params)  # You'll need to implement this method

        if self.response != None:
            response = self.response.function(**self.fit_result.params)
            best_fit = convolve1d(xs,response,0)
        residual = data_xs - best_fit

        # Plot styling
        color = kwargs.pop("color", "crimson")
        title = kwargs.pop("title", self.cross_section.name)
        ecolor = kwargs.pop("ecolor", "0.8")
        ms = kwargs.pop("ms", 2)

        # Top subplot: Data and fit
        ax[0].errorbar(wavelength, data_xs, 
                    # yerr=1./self.fit_result.weights,  # Assuming similar error handling
                    marker="o", 
                    color=color, 
                    ms=ms, 
                    zorder=-1, 
                    ecolor=ecolor, 
                    label="Cross-section data")
        
        ax[0].plot(wavelength, best_fit, color="0.4", label="Best fit")
        ax[0].plot(wavelength, xs, color="0.2", label="total xs")
        ax[0].set_ylabel("Total Cross Section [barn/sr]")
        ax[0].set_title(title)

        # Bottom subplot: Residuals
        ax[1].plot(wavelength, residual, color=color)
        ax[1].set_ylabel("Residuals [1σ]")
        ax[1].set_xlabel("λ [Å]")

        # Background plotting (if enabled)
        if plot_bg and self.background:
            self.background.plot(wl=wavelength, ax=ax[0], params=self.fit_result.params, **kwargs)
            ax[0].legend(["Cross-section data", "Background", "Total cross-section","Best fit"][::-1], 
                        fontsize=9, 
                        reverse=True, 
                        title=f"χ$^2$: {self.fit_result.redchi:.2f}")
        else:
            ax[0].legend(["Cross-section data","Total cross-section", "Best fit"][::-1], 
                        fontsize=9, 
                        reverse=True, 
                        title=f"χ$^2$: {self.fit_result.redchi:.2f}")

        # d-spacing plot (if enabled)
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
                        
                        # Label d-spacing lines
                        if len(self.cross_section.phases) > 1:
                            ax[0].text(dspace*2, dspace_label_pos, 
                                    f"{phase} {hkl}", 
                                    color="0.2", 
                                    zorder=-1, 
                                    fontsize=8, 
                                    transform=trans, 
                                    rotation=90, 
                                    va="top", 
                                    ha="right")
                        else:
                            ax[0].text(dspace*2, dspace_label_pos, 
                                    f"{hkl}", 
                                    color="0.2", 
                                    zorder=-1, 
                                    fontsize=8, 
                                    transform=trans, 
                                    rotation=90, 
                                    va="top", 
                                    ha="right")

        plt.subplots_adjust(hspace=0.05)
        return ax





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
        for i, material in enumerate(self._materials):
            mos = self._materials[material].get("mos", None)
            if mos: 
                param_name = f"η{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=mos, min=0.001, max=50, vary=vary)
                
                theta = self._materials[material].get("theta", 0.)
                param_name = f"θ{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=theta, min=0., max=180, vary=vary)
                
                phi = self._materials[material].get("phi", 0.)
                param_name = f"ϕ{i+1}"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=phi, min=0, max=360, vary=vary)
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
        for i, material in enumerate(self._materials):
            temp = self._materials[material].get("temp", None)
            if temp: 
                param_name = "temp"
                if param_name in self.params:
                    self.params[param_name].vary = vary
                else:
                    params.add(param_name, value=temp, min=77., max=1000, vary=vary)
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
        param_name = "thickness"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=1., min=0., max=5., vary=vary)
        
        param_name = "norm"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=1., min=0.1, max=10., vary=vary)
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
        param_name = "L0"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=L0, min=0.5, max=1.5, vary=vary)
        
        param_name = "t0"
        if param_name in self.params:
            self.params[param_name].vary = vary
        else:
            params.add(param_name, value=t0, min=-5e-6, max=5e6, vary=vary)
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
    
    def _make_lattice_params(self, vary: bool = False):
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

    
    def _make_extinction_params(self, vary: bool = False):
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

