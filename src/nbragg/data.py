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

        # Grouped data attributes
        self.is_grouped = False
        self.groups = None  # Dict mapping index -> table
        self.indices = None  # List of string indices
        self.group_shape = None  # Tuple (nx, ny) for 2D, (n,) for 1D, None for named

    def _normalize_index(self, index):
        """
        Normalize index for group lookup.
        Converts tuples like (10, 20) to strings like "(10,20)" for consistent access.
        Accepts both "(10,20)" and "(10, 20)" string formats.

        Parameters:
        -----------
        index : int, tuple, or str
            The index to normalize

        Returns:
        --------
        str
            String representation of the index (tuples without spaces)
        """
        if isinstance(index, tuple):
            # (10, 20) -> "(10,20)" (no spaces)
            return str(index).replace(" ", "")
        elif isinstance(index, str):
            # Remove spaces from string if it looks like a tuple: "(10, 20)" -> "(10,20)"
            return index.replace(" ", "")
        else:
            # 5 -> "5"
            return str(index)

    def _parse_string_index(self, string_idx):
        """
        Parse a string index back to its original form.
        "(10, 20)" -> (10, 20)
        "5" -> 5
        "center" -> "center"

        Parameters:
        -----------
        string_idx : str
            String representation of index

        Returns:
        --------
        tuple, int, or str
            Original index form
        """
        import ast
        try:
            # Try to parse as Python literal (for tuples and ints)
            parsed = ast.literal_eval(string_idx)
            return parsed
        except (ValueError, SyntaxError):
            # If parsing fails, it's a named string
            return string_idx

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
        
        # Ensure DataFrame has 'err' column with valid values
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
        elif df["err"].isna().all():
            # If err column exists but is all NaN (happens when reading 2-column CSV with 3 names)
            # calculate errors as sqrt of counts
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

    @classmethod
    def from_grouped(cls, signal, openbeam,
                     empty_signal: str = "", empty_openbeam: str = "",
                     tstep: float = 10.0e-6, L: float = 9,
                     L0: float = 1.0, t0: float = 0., dropna: bool = False,
                     pattern: str = "auto", indices: list = None, verbosity: int = 1,
                     n_jobs: int = -1):
        """
        Creates a Data object from grouped counts data using glob patterns.

        Supports 1D arrays, 2D grids, and named indices for spatially-resolved analysis.

        Parameters:
        -----------
        signal : str
            Glob pattern for signal files (e.g., "archive/pixel_*.csv" or "data/grid_*_x*_y*.csv").
            Can also be a folder path - all .csv files in the folder will be loaded.
        openbeam : str
            Glob pattern for openbeam files. Can also be a folder path.
        empty_signal : str, optional
            Glob pattern for empty signal files for background correction.
        empty_openbeam : str, optional
            Glob pattern for empty openbeam files for background correction.
        tstep : float, optional
            Time step (seconds) for converting time-of-flight to energy. Default is 10.0e-6.
        L : float, optional
            Distance (meters) used in the energy conversion from time-of-flight. Default is 9 m.
        L0 : float, optional
            Flight path scale factor from vary_tof optimization. Default is 1.0.
        t0 : float, optional
            Time offset correction (in tof units) from vary_tof optimization. Default is 0.
        dropna : bool, optional
            If True, remove rows with NaN values from data tables. Default is False.
        pattern : str, optional
            Coordinate extraction pattern. Default is "auto" which tries common patterns:
            - "x{x}_y{y}" for 2D grids (e.g., "grid_x10_y20.csv")
            - "idx{i}" or "pixel_{i}" for 1D arrays
            Custom patterns can use {x}, {y}, {i}, or {name}.
        indices : list, optional
            If provided, use these indices instead of extracting from filenames.
            Can be list of ints (1D), list of tuples (2D), or list of strings (named).
        verbosity : int, optional
            Verbosity level. If >= 1, shows progress bar. Default is 1.
        n_jobs : int, optional
            Number of parallel jobs for loading files. Default is -1 (use all CPUs).
            Set to 1 for sequential loading.

        Returns:
        --------
        Data
            A Data object with grouped data stored in self.groups.

        Examples:
        ---------
        # 2D grid from filenames like "pixel_x10_y20.csv"
        >>> data = Data.from_grouped("folder/pixel_*.csv", "folder_ob/pixel_*.csv")

        # 1D array with custom indices
        >>> data = Data.from_grouped("data/det_*.csv", "data_ob/det_*.csv", indices=[0, 1, 2, 3])

        # Named groups
        >>> data = Data.from_grouped("samples/*.csv", "ref/*.csv", indices=["sample1", "sample2"])
        """
        import glob
        import re
        import os

        # Find all matching files (support folder input)
        if os.path.isdir(signal):
            signal_files = sorted(glob.glob(os.path.join(signal, "*.csv")))
        else:
            signal_files = sorted(glob.glob(signal))

        if os.path.isdir(openbeam):
            openbeam_files = sorted(glob.glob(os.path.join(openbeam, "*.csv")))
        else:
            openbeam_files = sorted(glob.glob(openbeam))

        if not signal_files:
            raise ValueError(f"No files found matching pattern: {signal}")
        if not openbeam_files:
            raise ValueError(f"No files found matching pattern: {openbeam}")
        if len(signal_files) != len(openbeam_files):
            raise ValueError(f"Mismatch: {len(signal_files)} signal files vs {len(openbeam_files)} openbeam files")

        # Handle empty beam files if provided
        empty_signal_files = []
        empty_openbeam_files = []
        use_single_empty = False  # Flag for single empty file reuse

        if empty_signal and empty_openbeam:
            empty_signal_files = sorted(glob.glob(empty_signal))
            empty_openbeam_files = sorted(glob.glob(empty_openbeam))

            # Allow single empty file to be reused for all groups
            if len(empty_signal_files) == 1 and len(empty_openbeam_files) == 1:
                use_single_empty = True
            elif len(empty_signal_files) != len(signal_files) or len(empty_openbeam_files) != len(signal_files):
                raise ValueError(
                    f"Empty file count mismatch: {len(empty_signal_files)} empty signal, "
                    f"{len(empty_openbeam_files)} empty openbeam vs {len(signal_files)} signal files. "
                    f"Provide either 1 empty file (reused for all) or one per signal file."
                )

        # Extract or use provided indices
        if indices is not None:
            # Convert numpy arrays to list
            import numpy as np
            if isinstance(indices, np.ndarray):
                indices = indices.tolist()

            # User-provided indices
            if len(indices) != len(signal_files):
                raise ValueError(f"Number of indices ({len(indices)}) must match number of files ({len(signal_files)})")
            extracted_indices = indices
        else:
            # Auto-extract from filenames
            extracted_indices = cls._extract_indices_from_filenames(signal_files, pattern)

        # Determine group dimensionality and shape BEFORE converting to strings
        group_shape, is_2d, is_1d = cls._determine_group_shape(extracted_indices)

        # Convert all indices to strings for consistent access
        # For 2D: (10, 20) -> "(10,20)" (no spaces)
        # For 1D: 5 -> "5"
        # For named: "center" -> "center"
        string_indices = []
        for idx in extracted_indices:
            if isinstance(idx, tuple):
                # Convert tuple to string without spaces: "(10,20)"
                string_indices.append(str(idx).replace(" ", ""))
            elif isinstance(idx, str):
                string_indices.append(idx)  # "center"
            else:
                string_indices.append(str(idx))  # "5"

        extracted_indices = string_indices

        # Create Data object
        self_data = cls()
        self_data.is_grouped = True
        self_data.indices = extracted_indices
        self_data.group_shape = group_shape
        self_data.groups = {}
        self_data.L = L
        self_data.tstep = tstep

        # Helper function to load a single group
        def load_single_group(i, idx):
            """Load a single group's data files."""
            sig_file = signal_files[i]
            ob_file = openbeam_files[i]

            # Handle empty files - use single file if available, otherwise per-group
            if use_single_empty:
                es_file = empty_signal_files[0]
                eo_file = empty_openbeam_files[0]
            else:
                es_file = empty_signal_files[i] if empty_signal_files else ""
                eo_file = empty_openbeam_files[i] if empty_openbeam_files else ""

            # Create individual Data object for this group
            group_data = cls.from_counts(
                signal=sig_file,
                openbeam=ob_file,
                empty_signal=es_file,
                empty_openbeam=eo_file,
                tstep=tstep,
                L=L,
                L0=L0,
                t0=t0,
                dropna=dropna
            )

            return idx, group_data.table

        # Load groups in parallel or sequentially
        if n_jobs == 1:
            # Sequential loading with progress bar
            if verbosity >= 1:
                try:
                    from tqdm.auto import tqdm
                    iterator = tqdm(enumerate(extracted_indices), total=len(extracted_indices),
                                   desc=f"Loading {len(extracted_indices)} groups")
                except ImportError:
                    iterator = enumerate(extracted_indices)
            else:
                iterator = enumerate(extracted_indices)

            for i, idx in iterator:
                idx, table = load_single_group(i, idx)
                self_data.groups[idx] = table
        else:
            # Parallel loading
            from joblib import Parallel, delayed

            # Create progress bar if needed
            if verbosity >= 1:
                try:
                    from tqdm.auto import tqdm
                    pbar = tqdm(total=len(extracted_indices), desc=f"Loading {len(extracted_indices)} groups")
                except ImportError:
                    pbar = None
            else:
                pbar = None

            # Load groups in parallel
            results = Parallel(n_jobs=n_jobs)(
                delayed(load_single_group)(i, idx)
                for i, idx in enumerate(extracted_indices)
            )

            # Store results
            for idx, table in results:
                self_data.groups[idx] = table
                if pbar is not None:
                    pbar.update(1)

            if pbar is not None:
                pbar.close()

        # Set first group as default table for compatibility
        self_data.table = self_data.groups[extracted_indices[0]]

        return self_data

    @classmethod
    def _extract_indices_from_filenames(cls, filenames, pattern):
        """Extract indices from filenames based on pattern."""
        import re
        import os

        indices = []

        # Auto-detect pattern if needed
        if pattern == "auto":
            # Try common patterns
            test_name = os.path.basename(filenames[0])

            # Try 2D patterns - look for _x or _y to avoid matching dimension specs like 16x16
            match_2d = re.search(r'_x(\d+).*_y(\d+)', test_name, re.IGNORECASE)
            if not match_2d:
                # Try without underscores but with word boundaries
                match_2d = re.search(r'\bx(\d+)[_\s].*\by(\d+)', test_name, re.IGNORECASE)

            if match_2d:
                pattern = "_x{x}_y{y}"
            else:
                # Try 1D patterns - look for trailing numbers or with keywords
                match_1d = re.search(r'(?:idx|pixel|det)[_\s]*(\d+)', test_name, re.IGNORECASE)
                if not match_1d:
                    # Try just trailing number before extension
                    match_1d = re.search(r'_(\d+)\.', test_name)

                if match_1d:
                    pattern = "idx{i}"
                else:
                    # No pattern found - use filenames as indices (for named ROIs)
                    pattern = "{name}"

        # Extract based on pattern
        for fname in filenames:
            basename = os.path.basename(fname)

            if "{x}" in pattern and "{y}" in pattern:
                # 2D grid pattern - try multiple patterns
                # First try with underscores
                match = re.search(r'_x(\d+).*_y(\d+)', basename, re.IGNORECASE)
                if not match:
                    # Try with word boundaries
                    match = re.search(r'\bx(\d+)[_\s].*\by(\d+)', basename, re.IGNORECASE)
                if not match:
                    # Try simple pattern as last resort
                    match = re.search(r'(?<![\dx])x(\d+).*(?<![\dx])y(\d+)', basename, re.IGNORECASE)

                if match:
                    x, y = int(match.group(1)), int(match.group(2))
                    indices.append((x, y))
                else:
                    raise ValueError(f"Could not extract x,y coordinates from: {basename}. "
                                   f"Filename should contain _x<num> and _y<num> patterns.")

            elif "{i}" in pattern:
                # 1D array pattern - try multiple approaches
                # First try with keywords
                match = re.search(r'(?:idx|pixel|det)[_\s]*(\d+)', basename, re.IGNORECASE)
                if not match:
                    # Try trailing number before extension
                    match = re.search(r'_(\d+)\.', basename)
                if not match:
                    # Last resort - any number in the filename (rightmost)
                    matches = re.findall(r'(\d+)', basename)
                    if matches:
                        match = type('obj', (object,), {'group': lambda self, n: matches[-1]})()

                if match:
                    indices.append(int(match.group(1)))
                else:
                    raise ValueError(f"Could not extract index from: {basename}")

            elif "{name}" in pattern:
                # Named groups - use filename without extension
                name = os.path.splitext(basename)[0]
                indices.append(name)

            else:
                # Unknown pattern - use sequential
                indices.append(len(indices))

        return indices

    @classmethod
    def _determine_group_shape(cls, indices):
        """Determine group shape and dimensionality from indices."""
        # Check for empty indices (works with both lists and numpy arrays)
        if len(indices) == 0:
            return None, False, False

        first_idx = indices[0]

        # Check if 2D (tuples)
        if isinstance(first_idx, tuple) and len(first_idx) == 2:
            # 2D grid - use max coordinates + 1 to handle sparse grids
            xs = [idx[0] for idx in indices]
            ys = [idx[1] for idx in indices]
            return (max(ys) + 1, max(xs) + 1), True, False

        # Check if 1D (ints)
        elif isinstance(first_idx, (int, int.__class__)):
            # 1D array
            return (len(indices),), False, True

        # Named indices (strings)
        else:
            return None, False, False

    def __add__(self, other):
        """
        Adds two Data objects together by combining their signal and openbeam counts,
        then recalculating transmission with improved statistics.

        For grouped data, adds corresponding groups together to create a new grouped Data object.

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
            If grouped data don't have matching indices.
        TypeError
            If the objects don't have the necessary attributes for addition.
        """
        # Check that L and tstep are identical
        if self.L != other.L:
            raise ValueError(f"Cannot add Data objects with different L values: {self.L} != {other.L}")
        if self.tstep != other.tstep:
            raise ValueError(f"Cannot add Data objects with different tstep values: {self.tstep} != {other.tstep}")

        # Handle grouped data addition
        if self.is_grouped and other.is_grouped:
            # Check that indices match
            if set(self.indices) != set(other.indices):
                raise ValueError(
                    f"Cannot add grouped Data objects with different indices.\n"
                    f"self indices: {self.indices}\n"
                    f"other indices: {other.indices}"
                )

            # Create new grouped Data object
            result = Data()
            result.is_grouped = True
            result.indices = self.indices.copy()
            result.group_shape = self.group_shape
            result.groups = {}
            result.L = self.L
            result.tstep = self.tstep

            # Add each group pair
            for idx in self.indices:
                # Get tables for this index from both objects
                table1 = self.groups[idx]
                table2 = other.groups[idx]

                # Simple addition: combine transmission values
                # (Ideally we'd combine counts, but grouped data may not have signal/openbeam stored)
                combined_table = table1.copy()
                combined_table['trans'] = (table1['trans'] + table2['trans']) / 2  # Average
                # Error propagation for averaged values
                if 'err' in table1.columns and 'err' in table2.columns:
                    combined_table['err'] = np.sqrt(table1['err']**2 + table2['err']**2) / 2

                result.groups[idx] = combined_table

            return result

        # Handle non-grouped data addition
        elif not self.is_grouped and not other.is_grouped:
            # Validate that both objects have the necessary attributes
            if self.signal is None or self.openbeam is None:
                raise TypeError("Cannot add Data objects: this object was not created with from_counts()")
            if other.signal is None or other.openbeam is None:
                raise TypeError("Cannot add Data objects: other object was not created with from_counts()")

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

        else:
            # Mismatch between grouped and non-grouped
            raise TypeError("Cannot add grouped and non-grouped Data objects")

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

    def plot(self, index=None, **kwargs):
        """
        Plots the transmission data with error bars.

        Parameters:
        -----------
        index : int, tuple, or str, optional
            For grouped data, specify which group to plot:
            - int: 1D array index
            - tuple: (x, y) for 2D grid
            - str: named index
            If None and data is grouped, plots first group.
            If None and data is not grouped, plots the main table.
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

        # Determine which table to plot and set label
        plot_label = kwargs.pop("label", None)
        if self.is_grouped:
            # For grouped data
            if index is None:
                # Default to first group
                index = self.indices[0]
            # Normalize index for lookup (supports tuple, int, or string access)
            normalized_index = self._normalize_index(index)
            if normalized_index not in self.groups:
                raise ValueError(f"Index {index} not found in groups. Available indices: {self.indices}")
            table_to_plot = self.groups[normalized_index]
            # Add index to label if not provided
            if plot_label is None:
                plot_label = f"Index {index}"
        else:
            # For non-grouped data
            if index is not None:
                raise ValueError("Cannot specify index for non-grouped data")
            table_to_plot = self.table

        # Plot the data with error bars
        ax = table_to_plot.dropna().plot(x="wavelength",y="trans", yerr="err",
                                         xlim=xlim, ylim=ylim, logx=logx, ecolor=ecolor,
                                         xlabel=xlabel, ylabel=ylabel, label=plot_label, **kwargs)

        # Add legend if label was set
        if plot_label is not None:
            ax.legend()

        return ax

    def plot_map(self, wlmin=1.0, wlmax=5.0, **kwargs):
        """
        Plot transmission map averaged over wavelength range for grouped data.

        Parameters:
        -----------
        wlmin : float, optional
            Minimum wavelength for averaging (default: 1.0 Å).
        wlmax : float, optional
            Maximum wavelength for averaging (default: 5.0 Å).
        **kwargs : dict, optional
            Additional plotting parameters:
            - cmap : str, optional
              Colormap for 2D maps (default: 'viridis').
            - title : str, optional
              Plot title (default: auto-generated).
            - vmin, vmax : float, optional
              Color scale limits for 2D maps.
            - figsize : tuple, optional
              Figure size (width, height) in inches.

        Returns:
        --------
        matplotlib.Axes
            The axes of the plot.

        Raises:
        -------
        ValueError
            If called on non-grouped data.

        Examples:
        ---------
        >>> # For 2D grid data
        >>> data = Data.from_grouped("pixel_x*_y*.csv", "ob_x*_y*.csv", L=10, tstep=10e-6)
        >>> data.plot_map(wlmin=2.0, wlmax=4.0)

        >>> # For 1D array data
        >>> data = Data.from_grouped("pixel_*.csv", "ob_*.csv", L=10, tstep=10e-6)
        >>> data.plot_map(wlmin=1.0, wlmax=3.0)
        """
        import matplotlib.pyplot as plt
        import numpy as np

        if not self.is_grouped:
            raise ValueError("plot_map only works for grouped data")

        # Extract kwargs
        cmap = kwargs.pop("cmap", "viridis")
        title = kwargs.pop("title", None)
        vmin = kwargs.pop("vmin", None)
        vmax = kwargs.pop("vmax", None)
        figsize = kwargs.pop("figsize", None)

        # Calculate average transmission for each group
        avg_trans = {}
        for idx in self.indices:
            table = self.groups[idx]
            mask = (table['wavelength'] >= wlmin) & (table['wavelength'] <= wlmax)
            avg_trans[idx] = table.loc[mask, 'trans'].mean()

        # Create visualization based on group_shape
        if self.group_shape and len(self.group_shape) == 2:
            # 2D pcolormesh for proper block sizing
            # Extract unique x and y coordinates by parsing string indices
            xs = []
            ys = []
            for idx_str in self.indices:
                idx = self._parse_string_index(idx_str)
                if isinstance(idx, tuple) and len(idx) == 2:
                    xs.append(idx[0])
                    ys.append(idx[1])
            xs = sorted(set(xs))
            ys = sorted(set(ys))

            # Calculate grid spacing (block size)
            x_spacing = xs[1] - xs[0] if len(xs) > 1 else 1
            y_spacing = ys[1] - ys[0] if len(ys) > 1 else 1

            # Create coordinate arrays including edges for pcolormesh
            # Add half-spacing to create cell edges
            x_edges = np.array(xs) - x_spacing / 2
            x_edges = np.append(x_edges, xs[-1] + x_spacing / 2)
            y_edges = np.array(ys) - y_spacing / 2
            y_edges = np.append(y_edges, ys[-1] + y_spacing / 2)

            # Create 2D array for values
            trans_array = np.full((len(ys), len(xs)), np.nan)

            # Map indices to array positions
            x_map = {x: i for i, x in enumerate(xs)}
            y_map = {y: i for i, y in enumerate(ys)}

            for idx_str in self.indices:
                idx = self._parse_string_index(idx_str)
                if isinstance(idx, tuple) and len(idx) == 2:
                    x, y = idx
                    if x in x_map and y in y_map:
                        trans_array[y_map[y], x_map[x]] = avg_trans[idx_str]

            fig, ax = plt.subplots(figsize=figsize)
            im = ax.pcolormesh(x_edges, y_edges, trans_array, cmap=cmap, vmin=vmin, vmax=vmax,
                              shading='flat', **kwargs)
            ax.set_xlabel("X coordinate")
            ax.set_ylabel("Y coordinate")
            ax.set_aspect('equal')
            if title is None:
                title = f"Average Transmission Map ({wlmin:.1f}-{wlmax:.1f} Å)"
            ax.set_title(title)
            plt.colorbar(im, ax=ax, label="Transmission")
            return ax

        elif self.group_shape and len(self.group_shape) == 1:
            # 1D line plot - parse string indices back to integers
            indices_array = np.array([self._parse_string_index(idx) for idx in self.indices])
            trans_values = np.array([avg_trans[idx] for idx in self.indices])

            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(indices_array, trans_values, 'o-', **kwargs)
            ax.set_xlabel("Pixel index")
            ax.set_ylabel("Average Transmission")
            if title is None:
                title = f"Average Transmission ({wlmin:.1f}-{wlmax:.1f} Å)"
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            return ax

        else:
            # Bar chart for named indices
            fig, ax = plt.subplots(figsize=figsize)
            positions = np.arange(len(self.indices))
            trans_values = [avg_trans[idx] for idx in self.indices]

            ax.bar(positions, trans_values, **kwargs)
            ax.set_xticks(positions)
            ax.set_xticklabels(self.indices, rotation=45, ha='right')
            ax.set_ylabel("Average Transmission")
            if title is None:
                title = f"Average Transmission ({wlmin:.1f}-{wlmax:.1f} Å)"
            ax.set_title(title)
            plt.tight_layout()
            return ax
