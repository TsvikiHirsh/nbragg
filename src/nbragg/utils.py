import NCrystal as NC
from pathlib import Path
import numpy as np
import pickle
import json
import os

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

def get_cache_path():
    """
    Get the path to the materials cache file.
    Tries multiple locations in order of preference:
    1. User's home directory cache
    2. Package directory (if writable)
    3. Temporary directory as fallback
    """
    # Try user cache directory first (most reliable)
    try:
        if os.name == 'nt':  # Windows
            cache_dir = Path(os.environ.get('LOCALAPPDATA', Path.home() / 'AppData' / 'Local')) / 'nbragg'
        else:  # Linux/Mac
            cache_dir = Path.home() / '.cache' / 'nbragg'
        
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = cache_dir / 'materials_cache.pkl'
        
        # Test if writable
        cache_file.touch(exist_ok=True)
        return cache_file
    except (OSError, PermissionError):
        pass
    
    # Try package directory
    try:
        pkg_dir = Path(__file__).parent
        cache_file = pkg_dir / '.materials_cache.pkl'
        cache_file.touch(exist_ok=True)
        return cache_file
    except (OSError, PermissionError, NameError):
        pass
    
    # Fallback to temp directory
    import tempfile
    cache_dir = Path(tempfile.gettempdir()) / 'nbragg_cache'
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / 'materials_cache.pkl'

def get_user_materials_path():
    """
    Get the path to user-registered materials storage.
    This is separate from the cache and stores custom materials.
    """
    try:
        if os.name == 'nt':  # Windows
            data_dir = Path(os.environ.get('APPDATA', Path.home() / 'AppData' / 'Roaming')) / 'nbragg'
        else:  # Linux/Mac
            data_dir = Path.home() / '.local' / 'share' / 'nbragg'
        
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir / 'user_materials.json'
    except (OSError, PermissionError):
        # Fallback to temp
        import tempfile
        data_dir = Path(tempfile.gettempdir()) / 'nbragg_data'
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir / 'user_materials.json'

def save_cache(materials_dict):
    """Save materials dictionary to cache file."""
    try:
        cache_path = get_cache_path()
        with open(cache_path, 'wb') as f:
            pickle.dump(materials_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        return True
    except Exception as e:
        import warnings
        warnings.warn(f"Could not save materials cache: {e}")
        return False

def load_cache():
    """Load materials dictionary from cache file."""
    try:
        cache_path = get_cache_path()
        if cache_path.exists():
            with open(cache_path, 'rb') as f:
                return pickle.load(f)
    except Exception as e:
        import warnings
        warnings.warn(f"Could not load materials cache: {e}")
    return None

def save_user_materials(user_materials_dict):
    """Save user-registered materials to persistent storage."""
    try:
        user_path = get_user_materials_path()
        # Convert to JSON-serializable format
        json_dict = {}
        for key, value in user_materials_dict.items():
            json_dict[key] = {k: v for k, v in value.items() if v is not None}
        
        with open(user_path, 'w') as f:
            json.dump(json_dict, f, indent=2)
        return True
    except Exception as e:
        import warnings
        warnings.warn(f"Could not save user materials: {e}")
        return False

def load_user_materials():
    """Load user-registered materials from persistent storage."""
    try:
        user_path = get_user_materials_path()
        if user_path.exists():
            with open(user_path, 'r') as f:
                json_dict = json.load(f)
            
            # Convert back to proper format with None values
            materials_dict = {}
            for key, value in json_dict.items():
                materials_dict[key] = {
                    'mat': value.get('mat'),
                    'temp': value.get('temp', 300.0),
                    'mos': value.get('mos'),
                    'dir1': value.get('dir1'),
                    'dir2': value.get('dir2'),
                    'dirtol': value.get('dirtol'),
                    'theta': value.get('theta'),
                    'phi': value.get('phi'),
                    'a': value.get('a'),
                    'b': value.get('b'),
                    'c': value.get('c'),
                    'ext_method': value.get('ext_method'),
                    'ext_l': value.get('ext_l'),
                    'ext_Gg': value.get('ext_Gg'),
                    'ext_L': value.get('ext_L'),
                    'ext_tilt': value.get('ext_tilt'),
                    'weight': value.get('weight', 1.0),
                }
                # Add metadata fields
                for k, v in value.items():
                    if k.startswith('_'):
                        materials_dict[key][k] = v
            
            return materials_dict
    except Exception as e:
        import warnings
        warnings.warn(f"Could not load user materials: {e}")
    return {}

def make_materials_dict():
    """
    Populate the materials dictionary based on available ncmat files.
    Each entry uses the filename (without .ncmat) as key and contains
    a material specification dictionary compatible with CrossSection.
    """
    mat_dict = {}
    ncmat_files = NC.browseFiles()
    
    for mat in ncmat_files:
        if mat.name.endswith(".ncmat"):
            fullname = mat.name
            name_without_ext = fullname.replace(".ncmat", "")
            
            # Parse filename components
            name = ""
            formula = ""
            space_group = ""
            splitname = name_without_ext.split("_")
            
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
            
            # Create material specification compatible with CrossSection
            mat_dict[name_without_ext] = {
                'mat': fullname,
                'temp': 300.0,
                'mos': None,
                'dir1': None,
                'dir2': None,
                'dirtol': None,
                'theta': None,
                'phi': None,
                'a': None,
                'b': None,
                'c': None,
                'ext_method': None,
                'ext_l': None,
                'ext_Gg': None,
                'ext_L': None,
                'ext_tilt': None,
                'weight': 1.0,
                # Metadata for convenience
                '_name': name,
                '_formula': formula,
                '_space_group': space_group
            }
    
    return mat_dict

def initialize_materials():
    """
    Initialize the global materials dictionary.
    Uses cache if available, otherwise builds and caches it.
    Also loads user-registered materials.
    """
    global materials
    
    # Try to load from cache first
    cached_materials = load_cache()
    
    if cached_materials is not None:
        materials.update(cached_materials)
    else:
        # Build materials dictionary from NCrystal files
        materials.update(make_materials_dict())
        # Save to cache for next time
        save_cache(materials)
    
    # Load and merge user-registered materials
    user_materials = load_user_materials()
    if user_materials:
        materials.update(user_materials)

def rebuild_cache(save_to_cache=True):
    """
    Force rebuild of the materials dictionary from NCrystal files.
    
    Parameters:
    save_to_cache (bool): Whether to save the rebuilt cache (default True)
    
    Returns:
    dict: The rebuilt materials dictionary
    """
    global materials
    
    # Rebuild from NCrystal files
    new_materials = make_materials_dict()
    
    # Load user materials
    user_materials = load_user_materials()
    if user_materials:
        new_materials.update(user_materials)
    
    # Update global dictionary
    materials.clear()
    materials.update(new_materials)
    
    # Save to cache if requested
    if save_to_cache:
        save_cache(make_materials_dict())  # Cache only NCrystal materials
    
    return materials

def register_material(filename, material_name=None, save_persistent=True):
    """
    Register a material from a file and add it to the materials dictionary.
    
    Parameters:
    filename (str or Path): Path to the .ncmat file
    material_name (str, optional): Name to use for the material. If None, uses filename.
    save_persistent (bool): Whether to save to persistent user materials storage (default True)
    """
    filename = Path(filename)
    
    # Read and register with NCrystal
    with open(filename, "r") as fid:
        file_content = fid.read()
        NC.registerInMemoryFileData(filename.name, file_content)
    
    # Determine material name
    if material_name is None:
        material_name = filename.stem  # filename without extension
    
    # Parse filename components for metadata
    name = material_name
    formula = ""
    space_group = ""
    splitname = material_name.split("_")
    
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
    
    # Create material specification
    new_material = {
        'mat': filename.name,
        'temp': 300.0,
        'mos': None,
        'dir1': None,
        'dir2': None,
        'dirtol': None,
        'theta': None,
        'phi': None,
        'a': None,
        'b': None,
        'c': None,
        'ext_method': None,
        'ext_l': None,
        'ext_Gg': None,
        'ext_L': None,
        'ext_tilt': None,
        'weight': 1.0,
        '_name': name,
        '_formula': formula,
        '_space_group': space_group,
        '_custom': True  # Mark as user-registered
    }
    
    # Add to global materials
    materials[material_name] = new_material
    
    # Save to persistent storage if requested
    if save_persistent:
        user_materials = load_user_materials()
        user_materials[material_name] = new_material
        save_user_materials(user_materials)

def register_material_from_dict(material_dict, save_persistent=True):
    """
    Register one or more materials in the nbragg materials database from a dictionary.
    
    Parameters:
    material_dict (dict): A dictionary containing the material information. 
    save_persistent (bool): Whether to save to persistent user materials storage (default True)
    
    The dictionary can have one of two structures:
    
    1. Single material (flat structure):
       {
           'mat': 'Fe_sg225_Iron-gamma.ncmat',
           'temp': 300.0,
           'mos': None,
           'weight': 1.0,
           ...
       }
       In this case, a 'name' or '_name' key should be provided for the dictionary key.
    
    2. Multiple materials (nested structure):
       {
           'material1': {
               'mat': 'Fe_sg225_Iron-gamma.ncmat',
               'temp': 300.0,
               'mos': None,
               'weight': 1.0,
               ...
           },
           'material2': {
               'mat': 'ZrF4-beta_sg84.ncmat',
               'temp': 300.0,
               ...
           }
       }
    
    Returns:
    dict: The updated materials that were added
    """
    updated_materials = {}
    
    # Check if this is a nested dictionary (all values are dicts)
    if all(isinstance(v, dict) for v in material_dict.values()):
        # Nested structure - add all materials
        for material_name, material_info in material_dict.items():
            # Ensure it has required CrossSection format
            if 'mat' not in material_info:
                raise ValueError(f"Material '{material_name}' missing required 'mat' key")
            
            # Set defaults for missing keys
            full_material = {
                'mat': material_info['mat'],
                'temp': material_info.get('temp', 300.0),
                'mos': material_info.get('mos', None),
                'dir1': material_info.get('dir1', None),
                'dir2': material_info.get('dir2', None),
                'dirtol': material_info.get('dirtol', None),
                'theta': material_info.get('theta', None),
                'phi': material_info.get('phi', None),
                'a': material_info.get('a', None),
                'b': material_info.get('b', None),
                'c': material_info.get('c', None),
                'ext_method': material_info.get('ext_method', None),
                'ext_l': material_info.get('ext_l', None),
                'ext_Gg': material_info.get('ext_Gg', None),
                'ext_L': material_info.get('ext_L', None),
                'ext_tilt': material_info.get('ext_tilt', None),
                'weight': material_info.get('weight', 1.0),
                '_custom': True  # Mark as user-registered
            }
            
            # Copy over any metadata fields (starting with _)
            for key, value in material_info.items():
                if key.startswith('_'):
                    full_material[key] = value
            
            updated_materials[material_name] = full_material
    else:
        # Flat structure - single material
        if 'mat' not in material_dict:
            raise ValueError("Material dictionary missing required 'mat' key")
        
        # Get the name for this material
        material_name = material_dict.get('_name') or material_dict.get('name')
        if not material_name:
            # Try to extract from mat filename
            mat_file = material_dict['mat']
            material_name = mat_file.replace('.ncmat', '')
        
        # Set defaults for missing keys
        full_material = {
            'mat': material_dict['mat'],
            'temp': material_dict.get('temp', 300.0),
            'mos': material_dict.get('mos', None),
            'dir1': material_dict.get('dir1', None),
            'dir2': material_dict.get('dir2', None),
            'dirtol': material_dict.get('dirtol', None),
            'theta': material_dict.get('theta', None),
            'phi': material_dict.get('phi', None),
            'a': material_dict.get('a', None),
            'b': material_dict.get('b', None),
            'c': material_dict.get('c', None),
            'ext_method': material_dict.get('ext_method', None),
            'ext_l': material_dict.get('ext_l', None),
            'ext_Gg': material_dict.get('ext_Gg', None),
            'ext_L': material_dict.get('ext_L', None),
            'ext_tilt': material_dict.get('ext_tilt', None),
            'weight': material_dict.get('weight', 1.0),
            '_custom': True  # Mark as user-registered
        }
        
        # Copy over any metadata fields
        for key, value in material_dict.items():
            if key.startswith('_'):
                full_material[key] = value
        
        updated_materials[material_name] = full_material
    
    # Update global materials dictionary
    materials.update(updated_materials)
    
    # Save to persistent storage if requested
    if save_persistent:
        user_materials = load_user_materials()
        user_materials.update(updated_materials)
        save_user_materials(user_materials)
    
    return updated_materials

def clear_user_materials():
    """
    Clear all user-registered materials from persistent storage.
    Reloads materials from cache (NCrystal files only).
    """
    global materials
    
    # Remove user materials file
    try:
        user_path = get_user_materials_path()
        if user_path.exists():
            user_path.unlink()
    except Exception as e:
        import warnings
        warnings.warn(f"Could not clear user materials: {e}")
    
    # Reload from cache only
    materials.clear()
    cached_materials = load_cache()
    if cached_materials:
        materials.update(cached_materials)
    else:
        # Rebuild if no cache
        rebuild_cache()

# Initialize the materials dictionary at module import
# This is now lazy - only builds when first accessed
_initialized = False

def _ensure_initialized():
    """Ensure materials dictionary is initialized (lazy loading)."""
    global _initialized
    if not _initialized:
        initialize_materials()
        _initialized = True

# Override materials dict access to ensure initialization
class LazyMaterialsDict(dict):
    """Dictionary that initializes materials on first access."""
    
    def __getitem__(self, key):
        _ensure_initialized()
        return super().__getitem__(key)
    
    def __iter__(self):
        _ensure_initialized()
        return super().__iter__()
    
    def __len__(self):
        _ensure_initialized()
        return super().__len__()
    
    def keys(self):
        _ensure_initialized()
        return super().keys()
    
    def values(self):
        _ensure_initialized()
        return super().values()
    
    def items(self):
        _ensure_initialized()
        return super().items()
    
    def get(self, key, default=None):
        _ensure_initialized()
        return super().get(key, default)

# Replace materials dict with lazy version
materials = LazyMaterialsDict()