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
        A dataframe containing λ (Angstroms), transmission, and error values.
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
    
    @classmethod
    def _read_counts(cls, filename="run2_graphite_00000/graphite.csv"):
        """
        Reads the counts data from a CSV file and calculates errors if not provided.
        
        Parameters:
        -----------
        filename : str, optional
            The path to the CSV file containing time-of-flight (tof) and counts data. 
            Default is 'run2_graphite_00000/graphite.csv'.
        
        Returns:
        --------
        df : pandas.DataFrame
            A DataFrame containing columns: 'tof', 'counts', and 'err'. Errors are calculated 
            as the square root of counts if not provided in the file.
        """
        df = pd.read_csv(filename, names=["tof", "counts", "err"], header=None, skiprows=1)
        
        # If no error values provided, calculate as sqrt of counts
        if all(df["err"].isnull()):
            df["err"] = np.sqrt(df["counts"])
        
        # Store label from filename (without path and extension)
        df.attrs["label"] = filename.split("/")[-1].rstrip(".csv")
        
        return df
    
    @classmethod
    def from_counts(cls, signal: str, openbeam: str,
                    empty_signal: str = "", empty_openbeam: str = "",
                    tstep: float = 1.56255e-9, L: float = 10.59):
        """
        Creates a Data object from signal and open beam counts data, calculates transmission, 
        and converts tof to λ using energy-λ conversion.
        
        Parameters:
        -----------
        signal : str
            Path to the CSV file containing the signal data (tof, counts, err).
        openbeam : str
            Path to the CSV file containing the open beam data (tof, counts, err).
        empty_signal : str, optional
            Path to the CSV file containing the empty signal data for background correction. Default is an empty string.
        empty_openbeam : str, optional
            Path to the CSV file containing the empty open beam data for background correction. Default is an empty string.
        tstep : float, optional
            Time step (seconds) for converting time-of-flight (tof) to energy. Default is 1.56255e-9.
        L : float, optional
            Distance (meters) used in the energy conversion from time-of-flight. Default is 10.59 m.
        
        Returns:
        --------
        Data
            A Data object containing transmission and λ data.
        """
        # Read signal and open beam counts
        signal = cls._read_counts(signal)
        openbeam = cls._read_counts(openbeam)
        
        # Convert tof to energy using provided time step and distance
        signal["energy"] = utils.time2energy(signal["tof"] * tstep, L)
        
        # Convert energy to λ (Angstroms)
        signal["λ"] = signal["energy"].apply(NC.ekin2wl)
        
        # Calculate transmission and associated error
        transmission = signal["counts"] / openbeam["counts"]
        err = transmission * np.sqrt((signal["err"] / signal["counts"])**2 + 
                                     (openbeam["err"] / openbeam["counts"])**2)
        
        # If background (empty) data is provided, apply correction
        if empty_signal and empty_openbeam:
            empty_signal = cls._read_counts(empty_signal)
            empty_openbeam = cls._read_counts(empty_openbeam)
            
            transmission *= empty_openbeam["counts"] / empty_signal["counts"]
            err = transmission * np.sqrt(
                (signal["err"] / signal["counts"])**2 + 
                (openbeam["err"] / openbeam["counts"])**2 +
                (empty_signal["err"] / empty_signal["counts"])**2 + 
                (empty_openbeam["err"] / empty_openbeam["counts"])**2
            )
        
        # Construct a dataframe for λ, transmission, and error
        df = pd.DataFrame({
            "λ": signal["λ"],
            "trans": transmission,
            "err": err
        })
        
        # Set the label attribute from the signal file
        df.attrs["label"] = signal.attrs["label"]
        
        # Create and return the Data object
        self_data = cls()
        self_data.table = df#.set_index("λ")  # Set λ as the index
        self_data.tgrid = signal["tof"]
        
        return self_data
    
    @classmethod
    def from_transmission(cls, filename: str):
        """
        Creates a Data object directly from a transmission data file containing energy, transmission, and error values.
        Converts energy to λ and sets λ as the index.
        
        Parameters:
        -----------
        filename : str
            Path to the file containing the transmission data (energy, transmission, error) separated by whitespace.
        
        Returns:
        --------
        Data
            A Data object with the transmission data loaded into a dataframe.
        """
        df = pd.read_csv(filename, names=["energy", "trans", "err"], header=None, 
                         skiprows=0, delim_whitespace=True)
        
        # Convert energy to λ (Angstroms)
        df["λ"] = df["energy"].apply(NC.ekin2wl)
        
        # Create Data object and assign the dataframe with λ as index
        self_data = cls()
        self_data.table = df#.set_index("λ")
        
        return self_data
    
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
              Label for the x-axis (default: "λ [Å]").
            - ylabel : str, optional
              Label for the y-axis (default: "Transmission").
            - logx : bool, optional
              Whether to plot the x-axis on a logarithmic scale (default: False).
        
        Returns:
        --------
        matplotlib.Axes
            The axes of the plot containing the transmission data.
        """
        xlim = kwargs.pop("xlim", (0.5, 10))  # Default to λ range in Å
        ylim = kwargs.pop("ylim", (0., 1.))
        ecolor = kwargs.pop("ecolor", "0.8")
        xlabel = kwargs.pop("xlabel", "λ [Å]")
        ylabel = kwargs.pop("ylabel", "Transmission")
        logx = kwargs.pop("logx", False)  # Default is linear scale for λ
        
        # Plot the data with error bars
        return self.table.dropna().plot(x="λ",y="trans", yerr="err",
                                        xlim=xlim, ylim=ylim, logx=logx, ecolor=ecolor,
                                        xlabel=xlabel, ylabel=ylabel, **kwargs)
