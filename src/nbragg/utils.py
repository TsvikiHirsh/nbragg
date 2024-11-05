import NCrystal as NC
from pathlib import Path
import numpy as np

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

# Initialize materials as a global dictionary
materials = {}

def make_materials_dict():
    # populate the materials dictionary based on the available ncmat files in your system
    mat_dict = {}
    ncmat_files = NC.browseFiles()
    for mat in ncmat_files:
        if mat.name.endswith(".ncmat"):
            fullname = mat.name
            name = ""
            formula = ""
            space_group = ""
            splitname = fullname.replace(".ncmat", "").split("_")
            if len(splitname) == 3:
                formula, space_group, name = splitname
            elif len(splitname) == 1:
                formula = splitname[0]
                name = formula
            elif len(splitname) == 2:
                formula, space_group = splitname
                name = formula
            elif len(splitname) > 3:
                formula, space_group, *names = splitname
                name = "_".join(names)
            else:
                continue
            mat_dict[fullname] = {
                "name": name,
                "mat": fullname,
                "formula": formula,
                "space_group": space_group
            }
    return mat_dict

def initialize_materials():
    global materials
    materials.update(make_materials_dict())

def register_material(filename):
    filename = Path(filename)
    with open(filename, "r") as fid:
        NC.registerInMemoryFileData(filename.name, fid.read())
    initialize_materials()

import nbragg

def register_material_from_dict(material_dict):
    """
    Register one or more materials in the nbragg materials database from a dictionary.
    
    Parameters:
    material_dict (dict): A dictionary containing the material information. The dictionary can have one of two structures:
    
    1. Flat structure:
       {
           'name': 'Iron-gamma',
           'mat': 'Fe_sg225_Iron-gamma.ncmat',
           'formula': 'Fe',
           'space_group': 'sg225'
       }
    
    2. Nested structure:
       {
           'gamma': {
               'mat': 'ZrF4-beta_sg84.ncmat',
               'temp': 300.0,
               'mos': None,
               'k': None,
               'l': None,
               'weight': 0.5
           },
           'beta': {
               'mat': 'Fe_sg225_Iron-gamma.ncmat',
               'temp': 300.0,
               'mos': None,
               'k': None,
               'l': None,
               'weight': 0.5
           }
       }
    
    Returns:
    dict: The updated nbragg.materials dictionary.
    """
    updated_materials = {}
    
    if all(isinstance(v, dict) for v in material_dict.values()):
        # Nested dictionary structure
        for material_name, material_info in material_dict.items():
            updated_materials[material_name] = material_info
    else:
        # Flat dictionary structure
        updated_materials[material_dict['name']] = material_dict
    
    materials.update(updated_materials)
    return updated_materials

# Initialize the materials dictionary at module import
initialize_materials()
