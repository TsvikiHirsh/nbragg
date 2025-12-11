from nbragg import utils
import pandas as pd
import numpy as np
import NCrystal as NC

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
    signal : pandas.DataFrame or None
        The signal counts data (tof, counts, err).
    openbeam : pandas.DataFrame or None
        The open beam counts data (tof, counts, err).
    L : float or None
        Distance (meters) used in the energy conversion from time-of-flight.
    tstep : float or None
        Time step (seconds) for converting time-of-flight to energy.
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
        self.signal = None
        self.openbeam = None
        self.L = None
        self.tstep = None
    
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
                    tstep: float = 10.0e-6, L: float = 9,
                    L0: float = 1.0, t0: float = 0., dropna: bool = False):
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
        L0 : float, optional
            Flight path scale factor from vary_tof optimization. Default is 1.0.
            Values > 1.0 indicate a longer path, < 1.0 a shorter path.
        t0 : float, optional
            Time offset correction (in tof units) from vary_tof optimization. Default is 0.
        dropna : bool, optional
            If True, remove rows with NaN values from the data table. Default is False.

        Returns:
        --------
        Data
            A Data object containing transmission and wavelength data.
        """
        # Read signal and open beam counts
        signal = cls._read_counts(signal)
        openbeam = cls._read_counts(openbeam)

        # Apply L0 and t0 corrections using the same formula as in TransmissionModel._tof_correction
        # dtof = (1.0 - L0) * tof + t0, then corrected_tof = tof + dtof
        dtof = (1.0 - L0) * signal["tof"] + t0
        corrected_tof = signal["tof"] + dtof

        # Convert tof to energy using corrected time and nominal distance
        signal["energy"] = utils.time2energy(corrected_tof * tstep, L)

        # Convert energy to wavelength (Angstroms)
        signal["wavelength"] = signal["energy"].apply(NC.ekin2wl)

        # Calculate transmission and associated error
        transmission = signal["counts"] / openbeam["counts"]
        err = transmission * np.sqrt((signal["err"] / signal["counts"])**2 +
                                    (openbeam["err"] / openbeam["counts"])**2)

        # If background (empty) data is provided, apply correction
        # Check if empty_signal/empty_openbeam are non-empty (not empty string, not None, not empty DataFrame)
        has_empty_signal = (isinstance(empty_signal, pd.DataFrame) and not empty_signal.empty) or \
                          (isinstance(empty_signal, str) and empty_signal)
        has_empty_openbeam = (isinstance(empty_openbeam, pd.DataFrame) and not empty_openbeam.empty) or \
                            (isinstance(empty_openbeam, str) and empty_openbeam)

        if has_empty_signal and has_empty_openbeam:
            empty_signal = cls._read_counts(empty_signal)
            empty_openbeam = cls._read_counts(empty_openbeam)

            transmission *= empty_openbeam["counts"] / empty_signal["counts"]
            err = transmission * np.sqrt(
                (signal["err"] / signal["counts"])**2 +
                (openbeam["err"] / openbeam["counts"])**2 +
                (empty_signal["err"] / empty_signal["counts"])**2 +
                (empty_openbeam["err"] / empty_openbeam["counts"])**2
            )

        # Construct a dataframe for wavelength, transmission, and error
        df = pd.DataFrame({
            "wavelength": signal["wavelength"],
            "trans": transmission,
            "err": err
        })

        # Set the label attribute from the signal file
        df.attrs["label"] = signal.attrs["label"]

        # Drop NaN values if requested
        if dropna:
            df = df.dropna()

        # Create and return the Data object
        self_data = cls()
        self_data.table = df
        self_data.tgrid = signal["tof"]
        self_data.signal = signal
        self_data.openbeam = openbeam
        self_data.L = L
        self_data.tstep = tstep

        return self_data

    @classmethod
    def from_transmission(cls, input_data, index: str = "wavelength", dropna: bool = False):
        """
        Creates a Data object directly from transmission data containing wavelength/energy, transmission, and error values.

        Parameters:
        -----------
        input_data : str or pandas.DataFrame
            Path to a file containing the transmission data (wavelength/energy, transmission, error) separated by whitespace,
            or a pandas DataFrame with the transmission data.
        index : str, optional
            Name of the first column. If "energy", values will be converted to wavelength. Default is "wavelength".
        dropna : bool, optional
            If True, remove rows with NaN values from the data table. Default is False.

        Returns:
        --------
        Data
            A Data object with the transmission data loaded into a dataframe.
        """
        # Handle both file paths and DataFrames
        if isinstance(input_data, str):
            df = pd.read_csv(input_data, sep=r"\s+")
            df.columns = [index, "trans", "err"]
        elif isinstance(input_data, pd.DataFrame):
            df = input_data.copy()
            # Use the provided column names or assume they're already correct
            if len(df.columns) == 3:
                df.columns = [index, "trans", "err"]
            elif set(["trans", "err"]).issubset(df.columns):
                # Already has trans and err, just ensure index column is named correctly
                if index not in df.columns and len(df.columns) >= 3:
                    # Rename the first column that's not trans or err
                    for col in df.columns:
                        if col not in ["trans", "err"]:
                            df = df.rename(columns={col: index})
                            break
        else:
            raise TypeError("input_data must be a string (file path) or a pandas DataFrame")

        # Convert energy to wavelength if needed
        if index == "energy":
            df["wavelength"] = df["energy"].apply(NC.ekin2wl)

        # Drop NaN values if requested
        if dropna:
            df = df.dropna()

        # Create Data object and assign the dataframe
        self_data = cls()
        self_data.table = df

        return self_data
    
    def __add__(self, other):
        """
        Adds two Data objects together by combining their signal and openbeam counts,
        then recalculating transmission with improved statistics.

        Parameters:
        -----------
        other : Data
            Another Data object to add to this one.

        Returns:
        --------
        Data
            A new Data object with combined counts and recalculated transmission.

        Raises:
        -------
        ValueError
            If L or tstep parameters differ between the two Data objects.
        TypeError
            If the objects don't have the necessary attributes for addition.
        """
        # Validate that both objects have the necessary attributes
        if self.signal is None or self.openbeam is None:
            raise TypeError("Cannot add Data objects: this object was not created with from_counts()")
        if other.signal is None or other.openbeam is None:
            raise TypeError("Cannot add Data objects: other object was not created with from_counts()")

        # Check that L and tstep are identical
        if self.L != other.L:
            raise ValueError(f"Cannot add Data objects with different L values: {self.L} != {other.L}")
        if self.tstep != other.tstep:
            raise ValueError(f"Cannot add Data objects with different tstep values: {self.tstep} != {other.tstep}")

        # Add signal counts
        combined_signal = self.signal.copy()
        combined_signal["counts"] = self.signal["counts"] + other.signal["counts"]
        combined_signal["err"] = np.sqrt(self.signal["err"]**2 + other.signal["err"]**2)

        # Add openbeam counts
        combined_openbeam = self.openbeam.copy()
        combined_openbeam["counts"] = self.openbeam["counts"] + other.openbeam["counts"]
        combined_openbeam["err"] = np.sqrt(self.openbeam["err"]**2 + other.openbeam["err"]**2)

        # Calculate transmission and error with combined counts
        transmission = combined_signal["counts"] / combined_openbeam["counts"]
        err = transmission * np.sqrt(
            (combined_signal["err"] / combined_signal["counts"])**2 +
            (combined_openbeam["err"] / combined_openbeam["counts"])**2
        )

        # Create new dataframe with combined results
        df = pd.DataFrame({
            "wavelength": combined_signal["wavelength"],
            "trans": transmission,
            "err": err
        })

        # Combine labels
        label1 = self.table.attrs.get("label", "data1")
        label2 = other.table.attrs.get("label", "data2")
        df.attrs["label"] = f"{label1}+{label2}"

        # Create new Data object
        result = Data()
        result.table = df
        result.tgrid = self.tgrid
        result.signal = combined_signal
        result.openbeam = combined_openbeam
        result.L = self.L
        result.tstep = self.tstep

        return result

    def __iadd__(self, other):
        """
        In-place addition of another Data object to this one.

        Parameters:
        -----------
        other : Data
            Another Data object to add to this one.

        Returns:
        --------
        Data
            This Data object with updated values.

        Raises:
        -------
        ValueError
            If L or tstep parameters differ between the two Data objects.
        TypeError
            If the objects don't have the necessary attributes for addition.
        """
        result = self.__add__(other)
        self.table = result.table
        self.tgrid = result.tgrid
        self.signal = result.signal
        self.openbeam = result.openbeam
        self.L = result.L
        self.tstep = result.tstep
        return self

    def dropna(self, inplace=False):
        """
        Remove rows with NaN values from the data table.

        Parameters:
        -----------
        inplace : bool, optional
            If True, modify the Data object in place and return self.
            If False, return a new Data object with NaN values removed. Default is False.

        Returns:
        --------
        Data or None
            Returns a new Data object if inplace=False, or None if inplace=True.
        """
        if inplace:
            self.table = self.table.dropna()
            return self
        else:
            # Create a new Data object with dropna applied
            new_data = Data()
            new_data.table = self.table.dropna() if self.table is not None else None
            new_data.tgrid = self.tgrid
            new_data.signal = self.signal
            new_data.openbeam = self.openbeam
            new_data.L = self.L
            new_data.tstep = self.tstep
            return new_data

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
