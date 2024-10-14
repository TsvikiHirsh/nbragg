import numpy as np
import site
import os
from pathlib import Path
import shelve
from tqdm import tqdm
import requests
import mendeleev
import materials_compendium as mc
import pandas as pd

# Constants
SPEED_OF_LIGHT = 299792458  # m/s
MASS_OF_NEUTRON = 939.56542052 * 1e6 / (SPEED_OF_LIGHT ** 2)  # [eV s²/m²]


def time2energy(time, flight_path_length):
    """
    Convert time-of-flight to energy of the neutron.

    Parameters:
    time (float): Time-of-flight in seconds.
    flight_path_length (float): Flight path length in meters.

    Returns:
    float: Energy of the neutron in electronvolts (eV).
    """
    γ = 1 / np.sqrt(1 - (flight_path_length / time) ** 2 / SPEED_OF_LIGHT ** 2)
    return (γ - 1) * MASS_OF_NEUTRON * SPEED_OF_LIGHT ** 2  # eV

def energy2time(energy, flight_path_length):
    """
    Convert energy to time-of-flight of the neutron.

    Parameters:
    energy (float): Energy of the neutron in electronvolts (eV).
    flight_path_length (float): Flight path length in meters.

    Returns:
    float: Time-of-flight in seconds.
    """
    γ = 1 + energy / (MASS_OF_NEUTRON * SPEED_OF_LIGHT ** 2)
    return flight_path_length / SPEED_OF_LIGHT * np.sqrt(γ ** 2 / (γ ** 2 - 1))

def materials_dict():
    """
    Generate a dictionary containing information about materials, including their densities, 
    atomic densities, chemical formulas, and elemental/isotopic compositions.

    The function retrieves material data from the MaterialsCompendium package and organizes 
    it into a dictionary. Each material's total mass is calculated based on the atomic weights 
    and weight fractions of its constituent elements. The mass of each element reflects its 
    actual atomic mass, and the total mass is weighted by element fractions.

    Returns:
        materials (dict): A dictionary containing information for each material.
            Example structure:
            {
                'MaterialName': {
                    'name': 'MaterialName',
                    'density': <density>,  # g/cm³
                    'n': <atomic_density>,  # atoms per barn
                    'mass': <total_mass>,   # total mass (g/mol)
                    'formula': 'ChemicalFormula',
                    'elements': {
                        'ElementSymbol': {
                            'weight': <weight_fraction>,
                            'mass': <element_atomic_mass>,  # Actual atomic mass of the element (g/mol)
                            'isotopes': {
                                'IsotopeName': <isotope_weight_fraction>,
                                ...
                            }
                        }
                    }
                },
                ...
            }
    """
    materials = {}

    for material in mc.MaterialsCompendium:
        mat_name = material.Name
        density = material.Density
        n = material.MaterialAtomDensity
        formula = material.Formula
        total_mass = 0.0  # Initialize total weighted mass for the material

        elements = {}

        # Process each element in the material
        for element in material.Elements:
            name = element.Element
            weight_fraction = element.WeightFraction_whole  # Whole weight fraction of the element

            # Get the actual atomic mass of the element (e.g., 12 for Carbon)
            element_atomic_mass = element.AtomicMass  # Actual atomic mass (g/mol)

            # Add element to the elements dict with its actual atomic mass
            elements[name] = {
                "weight": weight_fraction,
                "mass": element_atomic_mass  # Actual atomic mass (g/mol)
            }

            # Calculate isotopic contributions if present
            isotopes = {}
            for isotope in element.Isotopes:
                iso_name = format_isotope(isotope.Isotope)
                iso_weight = isotope.WeightFraction
                if iso_weight > 0:
                    isotopes[iso_name] = iso_weight

            # Add isotopic data to the element
            elements[name]["isotopes"] = isotopes if isotopes else None  # If no isotopes, set None

            # Update total weighted mass for the material by summing element mass contributions weighted by fraction
            total_mass += element_atomic_mass * weight_fraction

        # Store the material information with calculated total weighted mass
        materials[mat_name] = {
            "name": mat_name,
            "density": density,   # g/cm³
            "n": n,               # atoms/barn
            "mass": total_mass,    # Total molar mass (g/mol), weighted by element fractions
            "formula": formula,
            "elements": elements   # Each element has its actual atomic mass (g/mol)
        }

    return materials





def elements_and_isotopes_dict():
    """
    Generate a dictionary of elements and isotopes with details such as 
    atoms per barn, total mass, and isotopic abundances.
    
    Returns:
        elements (dict): A dictionary containing information for each element.
        isotopes (dict): A dictionary containing information for each isotope.
    """
    elements = {}
    isotopes = {}
    
    import mendeleev
    all_elements = mendeleev.get_all_elements()

    for element in all_elements:
        name = element.name
        symbol = element.symbol

        # Handle cases where density or mass is missing (None)
        if element.density is None or element.mass is None:
            # Skip the element if either value is missing
            continue
        
        # Calculate atoms per barn and total mass
        atoms_per_barn = element.density / element.mass * 0.602214076  # Atoms/barn
        total_mass = element.mass  # Atomic mass (g/mol)
        
        # Initialize element entry with name, atoms per barn, and mass
        elements[name] = {
            "name": name,
            "n": atoms_per_barn,  # atoms per barn
            "mass": total_mass,   # total atomic mass (g/mol)
            "formula": symbol,    # Element symbol (formula)
            "density": element.density,  # g/cm³
            "elements": {symbol: {"weight": 1}}
        }

        # Add isotopic data if available
        element_isotopes = {f"{element.symbol}{iso.mass_number}": iso.abundance * 0.01
                            for iso in element.isotopes if iso.abundance is not None}
        
        if element_isotopes:
            elements[name]["elements"][symbol]["isotopes"] = element_isotopes
        
        # Add isotopes to a separate dictionary
        for iso in element.isotopes:
            if iso.abundance is not None:
                iso_name = f"{element.symbol}{iso.mass_number}"
                isotopes[iso_name] = {
                    "name": iso_name,
                    "n": element.density / iso.mass * 0.602214076,  # Atoms/barn for the isotope
                    "mass": iso.mass,  # Isotope mass
                    "formula": iso_name,
                    "density": element.density,
                    "elements": {symbol: {"weight": 1, "isotopes": {iso_name: 1.0}}}
                }

    return elements, isotopes




def format_isotope(isotope_string):
    # Find where the digits start in the string
    for i, char in enumerate(isotope_string):
        if char.isdigit():
            # Split the string into element and mass number parts
            element = isotope_string[:i]
            mass_number = isotope_string[i:]
            # Return the formatted string
            return f"{element}-{mass_number}"
    # Return the original string if no digits are found
    return isotope_string

def get_cache_path():
    # Get the user's site-packages directory
    user_site = site.getusersitepackages()
    # Create a subdirectory for our cache
    cache_dir = Path(user_site) / "nbragg"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir

def load_or_create_materials_cache():
    cache_path = get_cache_path() / "materials.db"
    
    if cache_path.exists():
        with shelve.open(str(cache_path.with_suffix(""))) as fid:
            materials = fid.get("materials")
            elements = fid.get("elements")
            isotopes = fid.get("isotopes")
        
        # If the cache is incomplete, regenerate it
        if materials is None or elements is None:
            return create_and_save_materials_cache()
        
        return materials, elements, isotopes
    else:
        return create_and_save_materials_cache()

def create_and_save_materials_cache():
    materials = materials_dict()
    elements, isotopes = elements_and_isotopes_dict()
    
    cache_path = get_cache_path() / "materials"
    with shelve.open(str(cache_path)) as fid:
        fid["materials"] = materials
        fid["elements"] = elements
        fid["isotopes"] = isotopes
    
    return materials, elements, isotopes

def download_xsdata():
    """Download the xsdata.npy file from GitHub and save it to the package directory."""
    url = 'https://github.com/lanl/trinidi-data/blob/main/xsdata.npy?raw=true'

    cache_path = get_cache_path() / "xsdata.npy"

    # Create the folder if it doesn't exist
    if not os.path.exists(cache_path.parent):
        os.makedirs(cache_path.parent)

    # Download with progress bar
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        total_size = int(response.headers.get('content-length', 0))
        with open(str(cache_path), 'wb') as f, tqdm(
            desc=f"Downloading {cache_path}",
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
            ncols=80
        ) as bar:
            for data in response.iter_content(chunk_size=1024):
                f.write(data)
                bar.update(len(data))
        print(f"File downloaded and saved to {cache_path}")
    else:
        raise Exception(f"Failed to download the file. Status code: {response.status_code}")
    

def interpolate(df):
    # Get underlying NumPy array and index values
    arr = df.values
    x = df.index.values  # Use the index values (energy bins) for interpolation

    # Loop through each column and interpolate NaNs using the index as x
    for j in range(arr.shape[1]):
        col = arr[:, j]
        valid_mask = ~np.isnan(col)
        
        if valid_mask.sum() > 1:  # Ensure there are at least two non-NaN values
            x_valid = x[valid_mask]     # Use index values for valid data points
            y_valid = col[valid_mask]   # Non-NaN values in the column
            x_nan = x[np.isnan(col)]    # Index values where NaNs are present
            
            if len(x_nan) > 0:  # Only interpolate if NaNs exist in the column
                # Interpolate NaNs using the actual index values
                col[np.isnan(col)] = np.interp(x_nan, x_valid, y_valid)

    # Return a DataFrame with the interpolated values
    return pd.DataFrame(arr, index=df.index, columns=df.columns)
