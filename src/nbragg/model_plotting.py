"""
Plotting methods for TransmissionModel, separated into a mixin class.
"""
from typing import TYPE_CHECKING, Optional, Union, List, Dict
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display

if TYPE_CHECKING:
    from nbragg.cross_section import CrossSection
    from nbragg.response import Response, Background
    import lmfit


class PlottingMixin:
    """
    Mixin class providing plotting functionality for TransmissionModel.

    This class contains all plotting-related methods that were originally
    part of the TransmissionModel class. It's designed to be mixed in with
    the main model class to provide visualization capabilities.

    Attributes (expected from parent class)
    ----------------------------------------
    cross_section : CrossSection
        The cross-section object containing phase data
    params : lmfit.Parameters
        Model parameters
    _materials : dict
        Dictionary of materials and their properties
    background : Background
        Background model (optional)
    response : Response
        Response function (optional)
    fit_result : lmfit.MinimizerResult
        Results from fitting (optional)
    fit_stages : list
        List of stage results from Rietveld fitting (optional)
    stages_summary : pandas.DataFrame
        Summary table of stage progression (optional)
    atomic_density : float
        Atomic density of the material

    Methods (expected from parent class)
    -------------------------------------
    eval(params, wl) : Calculate model at given wavelengths
    transmission(wl, **kwargs) : Calculate transmission
    """

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
