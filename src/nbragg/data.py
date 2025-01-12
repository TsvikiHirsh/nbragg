from nbragg import utils
import pandas as pd
import numpy as np
import NCrystal as NC
from scipy.interpolate import UnivariateSpline, RectBivariateSpline



class Data:
    """
    A class for handling neutron transmission data, including reading counts data, 
    calculating transmission, and plotting the results.
    
    Attributes:
    -----------
    table : pandas.DataFrame or None
        A dataframe containing wavelength (Angstroms), transmission, and error values.
    tgrid : pandas.Series or None
        A time-of-flight grid corresponding to the time steps in the data.
    """
    
    def __init__(self, **kwargs):
        """
        Initializes the Data object with optional keyword arguments.
        
        Parameters:
        -----------
        **kwargs : dict, optional
            Additional keyword arguments to set any instance-specific properties.
        """
        self.table = None
        self.tgrid = None
         self.original_table = None  # Store original data for rebinning
    
    @classmethod
    def _read_counts(cls, input_data, names=None):
        """
        Reads the counts data from a CSV file or a pandas DataFrame and calculates errors if not provided.
        
        Parameters:
        -----------
        input_data : str or pandas.DataFrame
            Either the path to the CSV file containing time-of-flight (tof) and counts data, 
            or a pandas DataFrame with the data.
        names : list, optional
            List of column names to use. If not provided, defaults to ["tof", "counts", "err"].
            Helps handle variations in column naming (e.g., "stacks" instead of "tof").
        
        Returns:
        --------
        df : pandas.DataFrame
            A DataFrame containing columns: 'tof', 'counts', and 'err'. Errors are calculated 
            as the square root of counts if not provided in the file.
        """
        # Default column names
        default_names = ["tof", "counts", "err"]
        
        # Process input based on type
        if isinstance(input_data, str):
            # If input is a file path, read CSV
            df = pd.read_csv(input_data, names=names or default_names, header=None, skiprows=1)
            # Store label from filename (without path and extension)
            df.attrs["label"] = input_data.split("/")[-1].rstrip(".csv")
        elif isinstance(input_data, pd.DataFrame):
            # If input is a DataFrame, create a copy to avoid modifying original
            df = input_data.copy()
            
            # Select the last 3 columns 
            if len(df.columns) > 3:
                df = df.iloc[:, -3:]
            
            # Process column names
            if names:
                # Rename columns if specific names are provided
                column_mapping = dict(zip(df.columns, names))
                df = df.rename(columns=column_mapping)
            else:
                # If no specific names, use default names 
                df.columns = default_names[:len(df.columns)]
            
            # Ensure we have the required columns
            required_columns = ["tof", "counts"]
            for col in required_columns:
                if col not in df.columns:
                    raise ValueError(f"DataFrame must contain a '{col}' column")
            
            # Use filename as label if available, otherwise use a default
            df.attrs["label"] = getattr(input_data, 'attrs', {}).get('label', 'input_data')
        else:
            raise TypeError("input_data must be a string (file path) or a pandas DataFrame")
        
        # Ensure DataFrame has 'err' column
        if "err" not in df.columns:
            # Try to find alternative error column names
            error_column_alternatives = ['error', 'std', 'std_dev', 'uncertainty']
            for alt_col in error_column_alternatives:
                if alt_col in df.columns:
                    df = df.rename(columns={alt_col: 'err'})
                    break
            else:
                # If no error column found, calculate as sqrt of counts
                df["err"] = np.sqrt(df["counts"])
        
        # Ensure consistent column order and names
        df = df[default_names[:len(df.columns)]]
        
        return df

    @classmethod
    def from_counts(cls, signal, openbeam,
                    empty_signal: str = "", empty_openbeam: str = "",
                    tstep: float = 10.0e-6, L: float = 9, sys_err: float = 0.):
        """
        Creates a Data object from signal and open beam counts data, calculates transmission, 
        and converts tof to wavelength using energy-wavelength conversion.
        
        Parameters:
        -----------
        signal : str or pandas.DataFrame
            Path to the CSV file or DataFrame containing the signal data (tof, counts, err).
        openbeam : str or pandas.DataFrame
            Path to the CSV file or DataFrame containing the open beam data (tof, counts, err).
        empty_signal : str or pandas.DataFrame, optional
            Path to the CSV file or DataFrame containing the empty signal data for background correction. 
            Default is an empty string.
        empty_openbeam : str or pandas.DataFrame, optional
            Path to the CSV file or DataFrame containing the empty open beam data for background correction. 
            Default is an empty string.
        tstep : float, optional
            Time step (seconds) for converting time-of-flight (tof) to energy. Default is 10.0e-6.
        L : float, optional
            Distance (meters) used in the energy conversion from time-of-flight. Default is 9 m.
        sys_err : float, optional
            Fractional systematical error to include in the transmission calculation 
            (e.g. sys_err=0.01 will include a 1% systematical error to the uncertainty)
        
        Returns:
        --------
        Data
            A Data object containing transmission and wavelength data.
        """
        # Read signal and open beam counts
        signal = cls._read_counts(signal)
        openbeam = cls._read_counts(openbeam)


        
        # Convert tof to energy using provided time step and distance
        signal["energy"] = utils.time2energy(signal["tof"] * tstep, L)
        
        # Convert energy to wavelength (Angstroms)
        signal["wavelength"] = signal["energy"].apply(NC.ekin2wl)
        
        # Calculate transmission and associated error
        transmission = signal["counts"] / openbeam["counts"]
        err = transmission * np.sqrt((signal["err"] / signal["counts"])**2 + 
                                    (openbeam["err"] / openbeam["counts"])**2 + sys_err**2)
        
        # If background (empty) data is provided, apply correction
        if empty_signal and empty_openbeam:
            empty_signal = cls._read_counts(empty_signal)
            empty_openbeam = cls._read_counts(empty_openbeam)
            
            transmission *= empty_openbeam["counts"] / empty_signal["counts"]
            err = transmission * np.sqrt(
                (signal["err"] / signal["counts"])**2 + 
                (openbeam["err"] / openbeam["counts"])**2 +
                (empty_signal["err"] / empty_signal["counts"])**2 + 
                (empty_openbeam["err"] / empty_openbeam["counts"])**2 +
                sys_err**2
            )
        
        # Construct a dataframe for wavelength, transmission, and error
        df = pd.DataFrame({
            "wavelength": signal["wavelength"],
            "trans": transmission,
            "err": err
        })
        
        # Set the label attribute from the signal file
        df.attrs["label"] = signal.attrs["label"]
        
        # Create and return the Data object
        self_data = cls()
        self_data.table = df
        self_data.tgrid = signal["tof"]
        self_data.L = L
        self_data.tstep = tstep
        
        return self_data

    @classmethod
    def from_transmission(cls, filename: str, index: str ="wavelength", L:float = 9.):
        """
        Creates a Data object directly from a transmission data file containing energy, transmission, and error values.
        Converts energy to wavelength and sets wavelength as the index.
        
        Parameters:
        -----------
        filename : str
            Path to the file containing the transmission data (energy, transmission, error) separated by whitespace.
        index : str 
            Optional energy to wavelegth convertion: If the index is "energy" it will be converted to wavelength
        L: float - flight path length [m]
        
        Returns:
        --------
        Data
            A Data object with the transmission data loaded into a dataframe.
        """
        df = pd.read_csv(filename, delim_whitespace=True)
        df.columns = [index, "trans", "err"]
        
        if index=="energy":
        # Convert energy to wavelength (Angstroms)
            df["wavelength"] = df["energy"].apply(NC.ekin2wl)

        
        # Create Data object and assign the dataframe with wavelength as index
        self_data = cls()
        self_data.table = df#.set_index("wavelength")
        λstep = df["wavelength"].iloc[1] - df["wavelength"].iloc[0]
        self_data.L = L
        self_data.tstep = λstep/3956.034*self_data.L
        
        return self_data

    def rebin(self, binwidth):
        """
        Rebin the data to new time bins.
        
        Parameters
        ----------
        binwidth : float
            New bin width in microseconds
        """
        if self.original_table is None:
            self.original_table = self.table.copy()
            
        # Convert binwidth from microseconds to seconds
        binwidth_s = binwidth * 1e-6
        
        # Create new time grid
        tmax = self.tgrid.max()
        new_tgrid = np.arange(0, tmax + binwidth_s, binwidth_s)
        
        # Prepare data for rebinning
        old_edges = np.r_[0, self.tgrid.values]
        new_edges = new_tgrid
        
        # Rebin the data
        rebinned_trans = rebin(old_edges, self.original_table['trans'].values, new_edges)
        rebinned_err = np.sqrt(rebin(old_edges, self.original_table['err'].values**2, new_edges))
        
        # Calculate wavelength for new time grid
        energy = utils.time2energy(new_tgrid[1:], self.L)
        wavelength = np.array([NC.ekin2wl(e) for e in energy])
        
        # Update table with rebinned data
        self.table = pd.DataFrame({
            'wavelength': wavelength,
            'trans': rebinned_trans,
            'err': rebinned_err
        })
        
        # Update time grid
        self.tgrid = pd.Series(new_tgrid[1:])
        self.tstep = binwidth_s
    
    def plot(self, **kwargs):
        """
        Plots the transmission data with error bars.
        
        Parameters:
        -----------
        **kwargs : dict, optional
            Additional plotting parameters:
            - xlim : tuple, optional
              Limits for the x-axis (default: (0.5, 10)).
            - ylim : tuple, optional
              Limits for the y-axis (default: (0., 1.)).
            - ecolor : str, optional
              Error bar color (default: "0.8").
            - xlabel : str, optional
              Label for the x-axis (default: "wavelength [Å]").
            - ylabel : str, optional
              Label for the y-axis (default: "Transmission").
            - logx : bool, optional
              Whether to plot the x-axis on a logarithmic scale (default: False).
        
        Returns:
        --------
        matplotlib.Axes
            The axes of the plot containing the transmission data.
        """
        xlim = kwargs.pop("xlim", (0.5, 10))  # Default to wavelength range in Å
        ylim = kwargs.pop("ylim", (0., 1.))
        ecolor = kwargs.pop("ecolor", "0.8")
        xlabel = kwargs.pop("xlabel", "wavelength [Å]")
        ylabel = kwargs.pop("ylabel", "Transmission")
        logx = kwargs.pop("logx", False)  # Default is linear scale for wavelength
        
        # Plot the data with error bars
        return self.table.dropna().plot(x="wavelength",y="trans", yerr="err",
                                        xlim=xlim, ylim=ylim, logx=logx, ecolor=ecolor,
                                        xlabel=xlabel, ylabel=ylabel, **kwargs)



class BoundedUnivariateSpline(UnivariateSpline):
    """1D spline that returns a constant for x outside the specified domain."""
    def __init__(self, x, y, fill_value=0.0, **kwargs):
        self.bnds = [x[0], x[-1]]
        self.fill_value = fill_value
        UnivariateSpline.__init__(self, x, y, **kwargs)
        
    def is_outside_domain(self, x):
        x = np.asarray(x)
        return np.logical_or(x < self.bnds[0], x > self.bnds[1])
        
    def __call__(self, x):
        outside = self.is_outside_domain(x)
        return np.where(outside, self.fill_value, UnivariateSpline.__call__(self, x))
        
    def integral(self, a, b):
        below_dx = np.max([0., self.bnds[0]-a])
        above_dx = np.max([0., b-self.bnds[1]])
        outside_contribution = (below_dx + above_dx) * self.fill_value
        a_f = np.max([a, self.bnds[0]])
        b_f = np.min([b, self.bnds[1]])
        if a_f >= b_f:
            return outside_contribution
        else:
            return outside_contribution + UnivariateSpline.integral(self, a_f, b_f)

def rebin(x1, y1, x2, interp_kind=3):
    """
    Rebin histogram values y1 from old bin edges x1 to new edges x2.
    
    Parameters
    ----------
    x1 : array_like
        Old bin edges
    y1 : array_like
        Old histogram values
    x2 : array_like
        New bin edges
    interp_kind : int or str, optional
        Interpolation type (3 for cubic spline or 'piecewise_constant')
        
    Returns
    -------
    array_like
        Rebinned histogram values
    """
    if interp_kind == 'piecewise_constant':
        return rebin_piecewise_constant(x1, y1, x2)
    else:
        return rebin_spline(x1, y1, x2, interp_kind=interp_kind)

def rebin_spline(x1, y1, x2, interp_kind):
    """Rebin using spline interpolation."""
    m = y1.size
    n = x2.size - 1
    
    x1_mid = (x1[:-1] + x1[1:]) / 2
    xx = np.hstack([x1[0], x1_mid, x1[-1]])
    yy = np.hstack([y1[0], y1, y1[-1]])
    
    spline = BoundedUnivariateSpline(xx, yy, s=0., k=interp_kind)
    areas1 = np.array([spline.integral(x1[i], x1[i+1]) for i in range(m)])
    
    x1_in_x2 = x1[np.logical_and(x1 > x2[0], x1 < x2[-1])]
    indices = np.searchsorted(x2, x1_in_x2)
    subbin_edges = np.insert(x2, indices, x1_in_x2)
    
    subbin_areas = np.array([spline.integral(subbin_edges[i], subbin_edges[i+1]) 
                            for i in range(subbin_edges.size-1)])
    
    subbin_mid = (subbin_edges[:-1] + subbin_edges[1:]) / 2
    sub2old = np.searchsorted(x1, subbin_mid) - 1
    sub2new = np.searchsorted(x2, subbin_mid) - 1
    
    y2 = np.zeros(n)
    for i in range(subbin_mid.size):
        if sub2old[i] == -1 or sub2old[i] == x1.size-1:
            continue
        y2[sub2new[i]] += (y1[sub2old[i]] * subbin_areas[i] / areas1[sub2old[i]])
    
    return y2