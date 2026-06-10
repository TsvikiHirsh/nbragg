import numpy as np
import pandas as pd
from scipy.stats import exponnorm
import matplotlib.pyplot as plt
from scipy.special import erfc, erfcx
from scipy.stats import uniform
import lmfit
import inspect
import warnings
import nbragg.utils as utils

class Response:
    def __init__(self, kind="jorgensen", vary: bool = False, wlstep: float = 0.1, tstep: float = 10e-6):
        """
        Initializes the Response object with specified parameters.

        Parameters
        ----------
        kind : str
            The type of response function to use. Options:

            - ``'jorgensen'`` — Jorgensen back-to-back exponentials with a
              Gaussian core, built on a fixed ``Δλ`` grid with step
              ``wlstep``. Cheap, but the *physical* convolution width depends
              on the user's wl sampling because ``convolve1d`` is sample-wise.
            - ``'jorgensen_inv'`` — same Jorgensen profile, but built on a
              Δλ grid whose step equals the data's actual ``dwl``. The
              physical width set by α/β/σ is then **invariant** under the
              user's choice of wl resolution. Use this whenever you fit/
              predict across different wl grids, or want α/β to describe the
              instrument rather than the grid.
            - ``'full_jorgensen'`` — TOF-domain Jorgensen with d-spacing
              dependence (α₀+α₁/d, β₀+β₁/d⁴, σ²=σ₀²+(σ₁d)²+(σ₂d²)²),
              GSAS-style µs units.
            - ``'square'`` — boxcar in time (width in µs).
            - ``'square_jorgensen'`` — Jorgensen convolved with a square.
            - ``'none'`` — no response (identity).
        vary : bool
            If True, the parameters can vary during fitting. Default False.
        wlstep : float
            Wavelength step size for the ``Δλ`` grid (Å), used by
            ``jorgensen`` and as a fallback for ``jorgensen_inv`` when
            ``data_dwl`` is not supplied. Default 0.1.
        tstep : float
            Time step size for ``tgrid`` in seconds, used by ``square`` and
            ``full_jorgensen``. Default 10e-6 (10 µs).
        """
        self.wlstep = wlstep
        self.tstep = tstep
        self.params = lmfit.Parameters()
        self.Δλ = np.arange(-20, 20, self.wlstep)
        self.tgrid = np.arange(-0.005, 0.005, tstep)  # Grid for time-based response (-5ms to 5ms with 10µs bins)
        self.kind = kind

        # Choose the response function
        if kind == "jorgensen":
            self.function = self.jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max=10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max=10000, vary=vary)

        elif kind == "square":
            self.function = self.square_response
            self.params = lmfit.Parameters()
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max=5000, vary=vary)

        elif kind == "square_jorgensen":
            self.function = self.square_jorgensen_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max=10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max=10000, vary=vary)
            self.params.add(f"width", value=tstep*1e6, min=tstep*1e6, max=5000, vary=vary)

        elif kind == "full_jorgensen":
            self.function = self.full_jorgensen_response
            self.params = lmfit.Parameters()
            # Units (user-facing): α, β are decay rates in 1/µs; σ is a
            # Gaussian width in µs. tgrid is still in seconds — the response
            # function converts internally. α0/β1/σ0/σ2 default to 0 and are
            # usually held fixed; σ1 is the dominant resolution term and
            # should be refined first (see staged_preset()). Defaults give a
            # near-delta response (FWHM ~ 1 tgrid sample = 10 µs) so the
            # convolution starts as essentially identity.
            self.params.add(f"α0", value=0.,    min=0.,    max=1000., vary=False)
            self.params.add(f"α1", value=0.5,   min=0.,    max=1000., vary=vary)
            self.params.add(f"β0", value=0.5,   min=0.,    max=1000., vary=vary)
            self.params.add(f"β1", value=0.0,   min=0.,    max=1000., vary=False)
            self.params.add(f"σ0", value=0.0,   min=0.0,   max=1000., vary=False)
            self.params.add(f"σ1", value=0.001, min=1e-6,  max=1000., vary=vary)
            self.params.add(f"σ2", value=0.0,   min=0.0,   max=1000., vary=False)

        elif kind == "jorgensen_inv":
            # Like 'jorgensen' but the profile is evaluated on a Δλ grid whose
            # step matches the data's actual wavelength step (passed in as
            # data_dwl by TransmissionModel.transmission). Result: the
            # physical resolution width set by α/β/σ is independent of how
            # finely the user samples wl.
            self.function = self.jorgensen_inv_response
            self.params = lmfit.Parameters()
            self.params.add(f"α0", value=3.67, min=0.001, max=10000, vary=vary)
            self.params.add(f"β0", value=3.06, min=0.001, max=10000, vary=vary)
            self.params.add(f"σ0", value=0.0,  min=0.0,   max=10.,   vary=False)
            self.params.add(f"σ1", value=0.0,  min=0.0,   max=10.,   vary=False)
            self.params.add(f"σ2", value=0.0,  min=0.0,   max=10.,   vary=False)

        elif kind == "none":
            self.function = self.empty_response
        else:
            raise NotImplementedError(f"Response kind '{kind}' is not supported. Use 'jorgensen', 'jorgensen_inv', 'square', 'square_jorgensen', 'full_jorgensen', or 'none'.")

    def staged_preset(self):
        """
        Recommended staged-fit ordering for this response kind.

        Returns a dict suitable for assigning to `TransmissionModel.stages`,
        merging with other stages as needed. For full_jorgensen the order
        follows Jorgensen-profile fitting best practice: refine σ1 first
        (dominant resolution term), then β0 (decay), then α1 (rise).
        """
        if self.kind == "full_jorgensen":
            return {
                "response_sig1": ["σ1"],
                "response_beta": ["σ1", "β0"],
                "response_alpha": ["σ1", "β0", "α1"],
            }
        if self.kind in ("jorgensen", "jorgensen_inv"):
            return {"response": ["α0", "β0"]}
        return {}

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

    def square_response(self, width=10, center: bool = True, **kwargs):
        """
        Returns a square response in time with a given width [µs].
        """
        width = width * 1e-6  # Convert to seconds
        loc = -0.5 * width if center else 0.
        tof_response = uniform.pdf(self.tgrid, loc=loc, scale=width)
        tof_response /= np.sum(tof_response)
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
            return np.nan_to_num(profile, 0)

    def jorgensen_inv_response(self, α0=3.67, β0=3.06, σ0=0.0, σ1=0.0, σ2=0.0,
                               data_dwl=None, data_wl=None, **kwargs):
        """
        Jorgensen peak profile, **invariant under the user's choice of wl
        sampling** (hence ``_inv``).

        Why this exists — the bug it fixes
        ----------------------------------
        The simple ``jorgensen_response`` builds its kernel on the fixed
        internal grid ``self.Δλ`` (step ``wlstep``, default 0.1 Å) and is then
        convolved with the cross-section via ``scipy.ndimage.convolve1d``,
        which is **unit-agnostic**: it applies the kernel sample-by-sample.
        That means the **physical** width of the smearing is
        ``N_kernel_samples × data_dwl`` — it scales linearly with how finely
        the user sampled ``wl``. Two equivalent measurements at different wl
        resolutions therefore produce different "instrument resolutions" out
        of the model, and fitted α/β values are tied to the grid rather than
        the instrument.

        How this version fixes it
        -------------------------
        ``jorgensen_inv`` evaluates the Jorgensen profile on a Δλ grid built
        with ``step = data_dwl`` — the median spacing of the user's wl array
        (passed in from ``TransmissionModel.transmission`` as ``data_wl``).
        Then ``N_kernel × data_dwl`` is the same physical width regardless of
        sampling. α, β, σ now describe the *instrument*, not the grid.

        Trade-offs
        ----------
        * **Use this when** you need α/β/σ to mean the same thing across
          datasets sampled differently, or when you fit data at one
          resolution and predict at another.
        * α/β values tuned against ``jorgensen`` are **not** transferable —
          ``jorgensen`` was effectively producing ``α·data_dwl/wlstep``-scaled
          kernels, so to get the same convolved transmission with
          ``jorgensen_inv`` you need to re-tune α/β to your instrument's
          actual physical resolution.
        * Convolution cost grows like ``1/data_dwl`` (more kernel samples for
          finer grids); ``n_half`` is capped at 20000 to keep this cheap on
          ultra-fine grids.

        Parameters (Å⁻¹ for α, β; Å for σ — d-spacing fixed at 1 Å)
        -----------------------------------------------------------
        α0, β0 : float
            Rise and decay rates at d = 1 Å (matches simple ``jorgensen``).
        σ0, σ1, σ2 : float
            Gaussian width terms; ``σ² = σ0² + (σ1·d)² + (σ2·d²)²`` with d = 1.
        data_dwl : float, optional
            Wavelength step of the data the response will be convolved with.
            Supplied automatically by ``TransmissionModel.transmission``. If
            omitted, falls back to ``data_wl`` (if provided) or
            ``self.wlstep`` (the simple-jorgensen behavior).
        data_wl : array-like, optional
            Data wavelength array; ``data_dwl`` is then taken as
            ``median(diff(data_wl))``.
        """
        # Resolve the target step
        if data_dwl is None and data_wl is not None and len(np.atleast_1d(data_wl)) > 1:
            data_dwl = float(np.median(np.diff(np.atleast_1d(data_wl))))
        if data_dwl is None or data_dwl <= 0:
            data_dwl = self.wlstep

        # Build a symmetric Δλ grid with exact step = data_dwl, covering
        # several decay constants and Gaussian widths on each side.
        d = 1.0
        α = α0 + 0.0
        β = β0 + 0.0
        σ = np.sqrt(σ0**2 + (σ1 * d)**2 + (σ2 * d * d)**2)
        decay = max(1.0 / max(α, 1e-6), 1.0 / max(β, 1e-6))
        margin = max(15.0 * decay, 8.0 * σ, 30.0 * data_dwl)
        # Cap the grid size to keep convolution cheap
        n_half = min(int(np.ceil(margin / data_dwl)), 20000)
        x = np.arange(-n_half, n_half + 1) * data_dwl

        sqrt2 = np.sqrt(2)
        σ_safe = max(σ, 1e-12)
        σ2_val = σ_safe * σ_safe
        scale = α * β / 2 / (α + β)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = α / 2 * (α * σ2_val + 2 * x)
            v = β / 2 * (β * σ2_val - 2 * x)
            y = (α * σ2_val + x) / (sqrt2 * σ_safe)
            z = (β * σ2_val - x) / (sqrt2 * σ_safe)

            ey = erfc(y); ez = erfc(z)
            term1 = np.exp(u) * ey
            term2 = np.exp(v) * ez
            term1[ey == 0] = 0
            term2[ez == 0] = 0

            profile = scale * (term1 + term2)
            profile = np.nan_to_num(profile, nan=0.0, posinf=0.0, neginf=0.0)
            total = profile.sum()
            if total > 0:
                profile = profile / total
            return profile

    def square_jorgensen_response(self, α0, β0, width, **kwargs):
        """
        Combines square and Jorgensen responses (not implemented in this update).
        """
        # Placeholder: Implement if needed
        raise NotImplementedError("square_jorgensen_response not implemented.")

    def full_jorgensen_response(self, wl=4.0, α0=0., α1=0.5, β0=0.5, β1=0., σ0=0., σ1=0.001, σ2=0., **kwargs):
        """
        Jorgensen TOF peak profile with d-spacing-dependent parameters.

        Implements eq. 3.40 of Jorgensen et al. (1978):
            h(Δ, σ, α, β) = αβ / (2(α+β)) · (exp(u) erfc(y) + exp(v) erfc(z))
        with
            α = α0 + α1 / d        [1/µs]
            β = β0 + β1 / d^4      [1/µs]
            σ² = σ0² + (σ1·d)² + (σ2·d²)²    [µs²]

        Parameter units are GSAS-style: α, β in 1/µs and σ in µs. The profile
        is evaluated on self.tgrid (seconds) — the function converts units
        internally so the user-facing values stay in the 0.001–1000 range.

        Uses the identity exp(u)·erfc(y) = exp(−x²/(2σ²))·erfcx(y) on the y≥0
        side to avoid exp(α²σ²/2) overflow when σ is large; the y<0 side is
        already bounded so it uses the direct form.

        Returns a 1-D normalized profile defined on self.tgrid.
        """
        d = wl / 2.0
        d2 = d * d
        d4 = d2 * d2

        # Convert user-facing units (1/µs, µs) to seconds-based internal units
        # so they match self.tgrid (seconds): 1/µs → ×1e6 to get 1/s; µs → ×1e-6 to get s.
        α = (α0 + α1 / d) * 1e6
        β = (β0 + β1 / d4) * 1e6
        σ2_val = (σ0**2 + (σ1 * d)**2 + (σ2 * d2)**2) * 1e-12
        σ2_val = max(σ2_val, 1e-30)
        σ = np.sqrt(σ2_val)

        if α + β <= 0:
            raise ValueError(f"α + β must be positive (got α={α}, β={β})")

        x = self.tgrid
        sqrt2 = np.sqrt(2)
        inv_2σ2 = 0.5 / σ2_val
        rad2sigma = 1.0 / (sqrt2 * σ)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            y = (α * σ2_val + x) * rad2sigma
            z = (β * σ2_val - x) * rad2sigma
            gauss = np.exp(-x * x * inv_2σ2)

            # term1 = exp(u)·erfc(y); use gauss·erfcx(y) for y>=0 (stable),
            # direct form for y<0 (where exp(u) is bounded ≤ exp(-α²σ²/2)·2)
            term1 = np.empty_like(x)
            m1 = y >= 0
            term1[m1] = gauss[m1] * erfcx(y[m1])
            u_neg = (α / 2) * (α * σ2_val + 2 * x[~m1])
            term1[~m1] = np.exp(u_neg) * erfc(y[~m1])

            term2 = np.empty_like(x)
            m2 = z >= 0
            term2[m2] = gauss[m2] * erfcx(z[m2])
            v_neg = (β / 2) * (β * σ2_val - 2 * x[~m2])
            term2[~m2] = np.exp(v_neg) * erfc(z[~m2])

            scale = α * β / (2 * (α + β))
            profile = scale * (term1 + term2)
            profile = np.nan_to_num(profile, nan=0.0, posinf=0.0, neginf=0.0)

            total = profile.sum()
            if total > 0:
                profile = profile / total
            return profile

    def _apply_asymmetry_correction(self, profile, xgrid, d, α0, β0, truncation_factor):
        """
        Applies asymmetry correction to the profile, mimicking GSAS computeAsymmetry.
        
        Parameters:
        -----------
        profile : array-like
            Input profile to correct.
        xgrid : array-like
            Grid of x values (time-of-flight offsets).
        d : float
            d-spacing.
        α0 : list/tuple
            [α0, α1] parameters.
        β0 : list/tuple
            [β0, β1] parameters.
        truncation_factor : float
            Truncation factor for asymmetry convolution.
        
        Returns:
        --------
        numpy.ndarray
            Asymmetry-corrected profile.
        """
        # Compute asymmetry parameters
        α = α0[0] + α0[1] / d
        β = β0[0] + β0[1] / d**4
        
        # Initialize output
        new_profile = np.zeros_like(profile)
        
        # Apply β asymmetry (first term)
        total_asymmetry = np.abs(β)
        if total_asymmetry > 0 and total_asymmetry < 1:
            direction = -1 if β < 0 else 1
            log_truncation = np.log(truncation_factor)
            truncation_angle = np.abs(-log_truncation / total_asymmetry)
            
            for i, x in enumerate(xgrid):
                function = profile[i]
                normalization = 1.0
                j = i + direction
                while 0 <= j < len(xgrid) and np.abs(xgrid[j] - x) < truncation_angle:
                    difference = np.abs(xgrid[j] - x)
                    exp_asymmetry = np.exp(-difference * total_asymmetry)
                    function += profile[j] * exp_asymmetry
                    normalization += exp_asymmetry
                    j += direction
                new_profile[i] = function / normalization if normalization > 0 else profile[i]
        
        # Apply α asymmetry (second term)
        total_asymmetry = np.abs(α)
        if total_asymmetry > 0 and total_asymmetry < 1:
            direction = -1 if α < 0 else 1
            log_truncation = np.log(truncation_factor)
            truncation_angle = np.abs(-log_truncation / total_asymmetry)
            
            for i, x in enumerate(xgrid):
                function = new_profile[i] if new_profile[i] != 0 else profile[i]
                normalization = 1.0
                j = i + direction
                while 0 <= j < len(xgrid) and np.abs(xgrid[j] - x) < truncation_angle:
                    difference = np.abs(xgrid[j] - x)
                    exp_asymmetry = np.exp(-difference * total_asymmetry)
                    function += new_profile[j] * exp_asymmetry
                    normalization += exp_asymmetry
                    j += direction
                new_profile[i] = function / normalization if normalization > 0 else new_profile[i]
        
        return new_profile
        
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
        if self.kind == "jorgensen" or self.kind == "square_jorgensen":
            xlabel = kwargs.pop("xlabel", "wavelength [Angstrom]")
            df = pd.Series(y, index=self.Δλ, name="Response")
            df.plot(ax=ax, xlabel=xlabel, **kwargs)
        elif self.kind == "jorgensen_inv":
            xlabel = kwargs.pop("xlabel", "Δλ [Angstrom]")
            x = (np.arange(len(y)) - len(y) // 2) * self.wlstep
            df = pd.Series(y, index=x, name="Response")
            df.plot(ax=ax, xlabel=xlabel, **kwargs)
        elif self.kind == "full_jorgensen":
            xlabel = kwargs.pop("xlabel", "tof [usec]")
            df = pd.Series(y, index=self.tgrid * 1e6, name="Response")
            df.plot(ax=ax, xlabel=xlabel, **kwargs)
        elif self.kind == "square":
            xlabel = kwargs.pop("xlabel", "tof [usec]")
            df = pd.Series(y, index=self.tgrid, name="Response")
            df.plot(ax=ax, xlabel=xlabel, **kwargs)
        else:
            raise ValueError("response plotting is only supported for 'jorgensen', 'jorgensen_inv', 'full_jorgensen', 'square', and 'square_jorgensen' kinds")


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
                self.params.add(f"bg{i}", value=0., min=-1e6, max= 1e6, vary=vary)

        elif kind == "polynomial5":
            self.function = self.polynomial5_background
            for i in range(5):
                self.params.add(f"bg{i}", value=0., min=-1e6, max= 1e6, vary=vary)
        elif kind == "sample_dependent":
            self.function = self.polynomial3_background
            self.params.add(f"k", value=1., min=1., max= 10., vary=vary)
            for i in range(3):
                self.params.add(f"bg{i}", value=0., min=-1e6, max= 1e6, vary=vary)


        elif kind == "constant":
            self.function = self.constant_background
            self.params.add('bg0', value=0.0, min=-1e6,max=1e6,vary=vary)
        elif kind == "none":
            self.function = self.empty_background
        else:
            raise NotImplementedError(f"Background kind '{kind}' is not supported. Use 'none', 'constant', 'polynomial3', 'sample_dependent', or 'polynomial5'.")

    def empty_background(self, wl, **kwargs):
        """
        Returns a zero background array.
        """
        return np.zeros_like(wl)

    def constant_background(self, wl, bg0=0., **kwargs):
        """
        Generates a constant background.

        Parameters:
        wl (np.ndarray): Wavelength values.
        b0 (float): Constant background value.
        """
        bg = np.full_like(wl, bg0)
        return np.where(bg>0,bg,0.)

    def polynomial3_background(self, wl, bg0=0., bg1=1., bg2=0., **kwargs):
        """
        Computwls a third-degree polynomial background.

        Parameters:
        wl (np.ndarray): Energy values.
        bg0 (float): Constant term.
        bg1 (float): Linear term.
        bg2 (float): Quadratic term.
        """
        bg = bg0 + bg1 * np.sqrt(wl) + bg2 / np.sqrt(wl)
        return np.where(bg>0,bg,0.)

    def polynomial5_background(self, wl, bg0=0., bg1=1., bg2=0., bg3=0., bg4=0., **kwargs):
        """
        Computes a fifth-degree polynomial background.

        Parameters:
        wl (np.ndarray): Wavelegth values.
        bg0 (float): Constant term.
        bg1 (float): Linear term.
        bg2 (float): Quadratic term.
        bg3 (float): Cubic term.
        bg4 (float): Quartic term.
        """
        bg = bg0 + bg1 * np.sqrt(wl) + bg2 / np.sqrt(wl) + bg3 * wl + bg4 * wl**2
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
