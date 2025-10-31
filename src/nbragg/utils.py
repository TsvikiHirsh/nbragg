import NCrystal as NC
from pathlib import Path
import numpy as np
import pickle
import json
import os
import warnings

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
_initialized = False
_initializing_materials = False  # Recursion guard for initialize_materials

def get_cache_path():
    """
    Get the path to the materials cache file.
    Tries multiple locations in order of preference:
    1. User's home directory cache
    2. Package directory (if writable)
    3. Temporary directory as fallback
    """
    try:
        if os.name == 'nt':  # Windows
            cache_dir = Path(os.environ.get('LOCALAPPDATA', Path.home() / 'AppData' / 'Local')) / 'nbragg'
        else:  # Linux/Mac
            cache_dir = Path(os.environ.get('HOME', '/tmp')) / '.cache' / 'nbragg'
        
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

def get_ncmat_storage_dir():
    """
    Get the directory for storing virtual NCMAT file contents.
    This allows materials to persist even if original files are deleted.
    """
    try:
        if os.name == 'nt':  # Windows
            data_dir = Path(os.environ.get('APPDATA', Path.home() / 'AppData' / 'Roaming')) / 'nbragg' / 'ncmat_files'
        else:  # Linux/Mac
            data_dir = Path.home() / '.local' / 'share' / 'nbragg' / 'ncmat_files'
        
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir
    except (OSError, PermissionError):
        # Fallback to temp
        import tempfile
        data_dir = Path(tempfile.gettempdir()) / 'nbragg_data' / 'ncmat_files'
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir

def save_cache(materials_dict):
    """Save materials dictionary to cache file."""
    try:
        cache_path = get_cache_path()
        with open(cache_path, 'wb') as f:
            pickle.dump(materials_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        return True
    except Exception as e:
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
        warnings.warn(f"Could not load materials cache: {e}")
    return None

def save_user_materials(user_materials_dict):
    """Save user-registered materials to persistent storage."""
    try:
        user_path = get_user_materials_path()
        # Convert to JSON-serializable format
        json_dict = {}
        for key, value in user_materials_dict.items():
            json_dict[key] = {k: v for k, v in value.items() if v is not None and k != '_ncmat_content'}
        
        with open(user_path, 'w') as f:
            json.dump(json_dict, f, indent=2)
        return True
    except Exception as e:
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
                    'ext_dist': value.get('ext_dist'),
                    'sans': value.get('sans'),
                    'weight': value.get('weight', 1.0),
                }
                # Add metadata fields
                for k, v in value.items():
                    if k.startswith('_'):
                        materials_dict[key][k] = v
                
                # Load virtual NCMAT content if it exists
                if '_ncmat_stored' in value and value['_ncmat_stored']:
                    ncmat_file = materials_dict[key]['mat']
                    stored_path = get_ncmat_storage_dir() / ncmat_file
                    if stored_path.exists():
                        with open(stored_path, 'r') as f:
                            ncmat_content = f.read()
                            # Register with NCrystal
                            NC.registerInMemoryFileData(ncmat_file, ncmat_content)
                            materials_dict[key]['_ncmat_content'] = ncmat_content
            
            return materials_dict
    except Exception as e:
        warnings.warn(f"Could not load user materials: {e}")
    return {}

def save_ncmat_content(filename, content):
    """
    Save NCMAT file content to persistent storage.
    
    Parameters:
    -----------
    filename : str
        Name of the NCMAT file
    content : str
        Content of the NCMAT file
    
    Returns:
    --------
    bool : True if successful, False otherwise
    """
    try:
        storage_dir = get_ncmat_storage_dir()
        filepath = storage_dir / filename
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    except Exception as e:
        warnings.warn(f"Could not save NCMAT content for {filename}: {e}")
        return False

def get_material_url(material_name):
    """
    Generate a URL to the NCrystal material documentation if available.
    
    Parameters:
    -----------
    material_name : str
        Name of the material file (with or without .ncmat extension)
    
    Returns:
    --------
    str or None : URL to the material documentation, or None if not available
    """
    if not material_name.endswith('.ncmat'):
        material_name += '.ncmat'
    
    base_url = "https://raw.githubusercontent.com/wiki/mctools/ncrystal/datalib/"
    return f"{base_url}{material_name}.inspect.pdf"

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
            
            # Generate URL for material documentation
            material_url = get_material_url(fullname)
            
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
                'ext_dist': None,
                'sans': None,
                'weight': 1.0,
                # Metadata for convenience
                '_name': name,
                '_formula': formula,
                '_space_group': space_group,
                '_url': material_url
            }
    
    return mat_dict


def initialize_materials():
    """Initialize the global materials dictionary, avoiding recursion issues."""
    global _initializing_materials, _initialized, materials

    # Recursion guard
    if _initializing_materials:
        raise RuntimeError(
            "Recursion detected in initialize_materials. "
            "Try calling initialize_materials() explicitly before accessing materials, "
            "or check for issues in get_cache_path() or IPython output handling."
        )

    _initializing_materials = True
    try:
        cached_materials = load_cache()

        # Always convert to a plain dict to avoid LazyMaterialsDict recursion
        if cached_materials is not None:
            base_data = dict(cached_materials)
        else:
            base_data = make_materials_dict()
            save_cache(base_data)

        # Bypass LazyMaterialsDict.update(), which calls _ensure_initialized()
        super(LazyMaterialsDict, materials).clear()
        super(LazyMaterialsDict, materials).update(base_data)

        # Mark initialization as complete before any access
        _initialized = True
        _initializing_materials = False

    finally:
        _initializing_materials = False


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

def register_material(filename=None, cif_source=None, save_persistent=True, store_content=True, **materials_dict):
    """
    Register one or more materials in the nbragg materials database.
    
    This function can handle multiple input formats:
    
    1. From a file:
       register_material("path/to/material.ncmat")
       register_material("path/to/material.ncmat", material_name="CustomName")
    
    2. From CIF file or COD ID:
       register_material(cif_source="path/to/material.cif", material_name="MyCIF")
       register_material(cif_source="codid::7123352", material_name="Bismuth")
    
    3. From a dictionary (single material):
       register_material(my_material={'mat': 'Fe.ncmat', 'temp': 500.0})
    
    4. From a dictionary (multiple materials):
       register_material(
           iron={'mat': 'Fe.ncmat', 'temp': 300.0},
           copper={'mat': 'Cu.ncmat', 'temp': 400.0}
       )
    
    5. From CrossSection.materials:
       xs = CrossSection(...)
       register_material(**xs.materials)
       # or with custom names:
       register_material(custom=xs.materials['phase1'])
    
    Parameters:
    -----------
    filename : str or Path, optional
        Path to an .ncmat file to register. If provided, other sources are ignored.
    cif_source : str, optional
        Either a path to a CIF file or a COD ID in format "codid::XXXXXX"
    save_persistent : bool, optional
        Whether to save to persistent user materials storage (default True)
    store_content : bool, optional
        Whether to store the NCMAT file content for offline access (default True)
        This ensures materials persist even if original files are deleted
    **materials_dict : dict
        Keyword arguments where each key is a material name and value is either:
        - A material specification dict with 'mat' key
        - A nested dict of material specifications (from CrossSection.materials)
    
    Returns:
    --------
    dict : The materials that were registered
    """
    updated_materials = {}
    
    # Handle CIF input
    if cif_source is not None:
        try:
            # Create NCMATComposer from CIF
            composer = NC.NCMATComposer(cif_source)
            
            # Determine material name
            material_name = materials_dict.pop('material_name', None)
            if material_name is None:
                if cif_source.startswith('codid::'):
                    material_name = f"COD_{cif_source.split('::')[1]}"
                else:
                    material_name = Path(cif_source).stem
            
            # Generate NCMAT filename
            ncmat_filename = f"{material_name}.ncmat"
            
            # Get NCMAT content
            ncmat_content = composer.create_ncmat()
            
            # Register with NCrystal
            NC.registerInMemoryFileData(ncmat_filename, ncmat_content)
            
            # Store content if requested
            if store_content:
                save_ncmat_content(ncmat_filename, ncmat_content)
            
            # Determine URL (if from COD)
            material_url = None
            if cif_source.startswith('codid::'):
                cod_id = cif_source.split('::')[1]
                material_url = f"http://www.crystallography.net/cod/{cod_id}.html"
            elif Path(cif_source).suffix == '.cif':
                material_url = f"file://{Path(cif_source).absolute()}"
            
            # Create material specification
            new_material = {
                'mat': ncmat_filename,
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
                'ext_dist': None,
                'sans': None,
                'weight': 1.0,
                '_name': material_name,
                '_formula': '',
                '_space_group': '',
                '_custom': True,
                '_cif_source': cif_source,
                '_url': material_url,
                '_ncmat_stored': store_content
            }
            
            if store_content:
                new_material['_ncmat_content'] = ncmat_content
            
            updated_materials[material_name] = new_material
            
        except Exception as e:
            raise RuntimeError(f"Failed to process CIF source '{cif_source}': {e}")
    
    # Handle file input
    elif filename is not None:
        filename = Path(filename)
        
        # Read file content
        with open(filename, "r") as fid:
            file_content = fid.read()
        
        # Register with NCrystal
        NC.registerInMemoryFileData(filename.name, file_content)
        
        # Store content if requested
        if store_content:
            save_ncmat_content(filename.name, file_content)
        
        # Determine material name
        material_name = materials_dict.pop('material_name', None)
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
            'ext_dist': None,
            'sans': None,
            'weight': 1.0,
            '_name': name,
            '_formula': formula,
            '_space_group': space_group,
            '_custom': True,
            '_url': f"file://{filename.absolute()}",
            '_ncmat_stored': store_content
        }
        
        if store_content:
            new_material['_ncmat_content'] = file_content
        
        updated_materials[material_name] = new_material
        
    
    # Handle dictionary input
    elif materials_dict:
        for material_name, material_info in materials_dict.items():
            # Check if material_info is a nested dictionary (like from CrossSection.materials)
            if isinstance(material_info, dict):
                # Check if this is a nested structure with multiple materials
                if all(isinstance(v, dict) and 'mat' in v for v in material_info.values()):
                    # This is a nested structure like {'phase1': {...}, 'phase2': {...}}
                    for sub_name, sub_info in material_info.items():
                        full_name = f"{material_name}_{sub_name}"
                        processed_mat = _process_material_dict(sub_info)
                        
                        # Try to store content if it's a custom material
                        if store_content and sub_info.get('mat'):
                            _try_store_material_content(processed_mat)
                        
                        updated_materials[full_name] = processed_mat
                else:
                    # This is a single material specification
                    if 'mat' not in material_info:
                        raise ValueError(f"Material '{material_name}' missing required 'mat' key")
                    
                    processed_mat = _process_material_dict(material_info)
                    
                    # Try to store content if it's a custom material
                    if store_content and material_info.get('mat'):
                        _try_store_material_content(processed_mat)
                    
                    updated_materials[material_name] = processed_mat
            else:
                raise ValueError(f"Material specification for '{material_name}' must be a dictionary")
    else:
        raise ValueError("Either 'filename', 'cif_source', or material dictionaries must be provided")
    
    # Update global materials dictionary
    materials.update(updated_materials)
    
    # Save to persistent storage if requested
    if save_persistent:
        user_materials = load_user_materials()
        user_materials.update(updated_materials)
        save_user_materials(user_materials)
    
    return {k: MaterialSpec(v) if not isinstance(v, MaterialSpec) else v 
        for k, v in updated_materials.items()}

def _try_store_material_content(material_dict):
    """
    Try to store the NCMAT content for a material if accessible.
    
    Parameters:
    -----------
    material_dict : dict
        Material dictionary to update with stored content info
    """
    try:
        mat_file = material_dict.get('mat')
        if mat_file:
            # Try to get the content from NCrystal
            text_data = NC.createTextData(mat_file)
            ncmat_content = text_data.rawData
            
            # Save the content
            if save_ncmat_content(mat_file, ncmat_content):
                material_dict['_ncmat_stored'] = True
                material_dict['_ncmat_content'] = ncmat_content
    except Exception:
        # If we can't access the content, that's okay
        pass

def _process_material_dict(material_info):
    """
    Process a material info dictionary into the standard format.
    
    Parameters:
    -----------
    material_info : dict
        Material specification dictionary
    
    Returns:
    --------
    dict : Processed material dictionary with all required keys
    """
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
        'ext_dist': material_info.get('ext_dist', None),
        'sans': material_info.get('sans', None),
        'weight': material_info.get('weight', 1.0),
        '_custom': True  # Mark as user-registered
    }
    
    # Copy over any metadata fields (starting with _)
    for key, value in material_info.items():
        if key.startswith('_') and key != '_custom':
            full_material[key] = value
    
    return full_material

def get_material_info(material_key, show_url=False):
    """
    Get information about a material in a formatted way.
    
    Parameters:
    -----------
    material_key : str
        Key of the material in the materials dictionary
    show_url : bool, optional
        Whether to display clickable URL (default False)
    
    Returns:
    --------
    str : Formatted material information
    """
    _ensure_initialized()
    
    if material_key not in materials:
        return f"Material '{material_key}' not found"
    
    mat = materials[material_key]
    
    info_lines = [f"Material: {material_key}"]
    info_lines.append(f"  File: {mat['mat']}")
    
    if mat.get('_name'):
        info_lines.append(f"  Name: {mat['_name']}")
    if mat.get('_formula'):
        info_lines.append(f"  Formula: {mat['_formula']}")
    if mat.get('_space_group'):
        info_lines.append(f"  Space Group: {mat['_space_group']}")
    
    info_lines.append(f"  Temperature: {mat['temp']} K")
    
    if mat.get('_cif_source'):
        info_lines.append(f"  CIF Source: {mat['_cif_source']}")
    
    if mat.get('_custom'):
        info_lines.append("  Status: User-registered")
    
    if mat.get('_ncmat_stored'):
        info_lines.append("  Storage: Content stored locally")
    
    if mat.get('_url'):
        if show_url:
            info_lines.append(f"  URL: {mat['_url']}")
        else:
            info_lines.append(f"  Documentation: Available (use show_url=True)")
    
    return '\n'.join(info_lines)

def list_materials(filter_str=None, custom_only=False):
    """
    List all available materials with optional filtering.
    
    Parameters:
    -----------
    filter_str : str, optional
        Filter materials by name/formula (case-insensitive)
    custom_only : bool, optional
        Only show user-registered materials (default False)
    
    Returns:
    --------
    list : List of material keys matching the criteria
    """
    _ensure_initialized()
    
    result = []
    for key, mat in materials.items():
        # Apply filters
        if custom_only and not mat.get('_custom'):
            continue
        
        if filter_str:
            search_str = filter_str.lower()
            if not any(search_str in str(mat.get(field, '')).lower() 
                      for field in ['_name', '_formula', 'mat', '_space_group']):
                continue
        
        result.append(key)
    
    return sorted(result)

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
        warnings.warn(f"Could not clear user materials: {e}")
    
    # Reload from cache only
    materials.clear()
    cached_materials = load_cache()
    if cached_materials:
        materials.update(cached_materials)
    else:
        # Rebuild if no cache
        rebuild_cache()

def _ensure_initialized():
    """Ensure materials dictionary is initialized (lazy loading)."""
    global _initialized
    if not _initialized:
        initialize_materials()
        _initialized = True

class Hyperlink:
    """
    A simple hyperlink representation that displays as 'hyperlink' in repr
    but can be clicked in Jupyter/IPython environments.
    """
    
    def __init__(self, url):
        self.url = url
    
    def __repr__(self):
        return self.url
    
    def _repr_html_(self):
        """
        Jupyter/IPython HTML representation - creates a clickable link.
        """
        return f'<a href="{self.url}" target="_blank">hyperlink</a>'
    
    def __str__(self):
        return self.url


class MaterialSpec(dict):
    """
    A dictionary subclass for material specifications that supports weight multiplication
    and provides enhanced representation with hyperlinks.
    
    Usage:
        material = nbragg.materials["Al2O3_sg167_Corundum"]
        weighted_material = material * 0.5  # Creates new MaterialSpec with weight=0.5
        xs = CrossSection(phase1=material * 0.3, phase2=other_material * 0.7)
    """
    
    def __mul__(self, weight):
        """
        Multiply material by a weight factor.
        
        Parameters:
        -----------
        weight : float
            Weight factor to apply to the material
        
        Returns:
        --------
        MaterialSpec : New MaterialSpec with updated weight
        """
        if not isinstance(weight, (int, float)):
            return NotImplemented
        
        # Create a copy of the material
        new_material = MaterialSpec(self)
        new_material['weight'] = float(weight)
        return new_material
    
    def __rmul__(self, weight):
        """
        Right multiplication (weight * material).
        
        Parameters:
        -----------
        weight : float
            Weight factor to apply to the material
        
        Returns:
        --------
        MaterialSpec : New MaterialSpec with updated weight
        """
        return self.__mul__(weight)
    
    def __repr__(self):
        """
        Enhanced representation with hyperlink for URL.
        """
        # Create a copy of the dict for display
        display_dict = dict(self)
        
        # Replace URL with hyperlink object if present
        if '_url' in display_dict and display_dict['_url']:
            url = display_dict['_url']
            display_dict['_url'] = Hyperlink(url)
        
        return repr(display_dict)
    
    def _repr_pretty_(self, p, cycle):
        """
        IPython/Jupyter pretty printing support.
        """
        if cycle:
            p.text('{...}')
            return
        
        with p.group(1, '{', '}'):
            if self:
                for idx, (key, value) in enumerate(self.items()):
                    if idx:
                        p.text(',')
                        p.breakable()
                    
                    p.text(repr(key))
                    p.text(': ')
                    
                    # Special handling for URL
                    if key == '_url' and value:
                        p.text(repr(Hyperlink(value)))
                    else:
                        p.pretty(value)


class LazyMaterialsDict(dict):
    """Dictionary that initializes materials on first access and returns MaterialSpec objects."""
    
    def __getitem__(self, key):
        _ensure_initialized()

        # Try exact match first
        try:
            material_dict = super().__getitem__(key)
        except KeyError:
            # If key ends with .ncmat, try without extension
            if key.endswith('.ncmat'):
                try:
                    material_dict = super().__getitem__(key[:-6])
                except KeyError:
                    raise KeyError(key)
            # If key doesn't end with .ncmat, try with extension
            else:
                try:
                    material_dict = super().__getitem__(f"{key}.ncmat")
                except KeyError:
                    raise KeyError(key)

        # Wrap in MaterialSpec if not already wrapped
        if not isinstance(material_dict, MaterialSpec):
            material_dict = MaterialSpec(material_dict)
            # Update the stored value so we don't re-wrap every time
            # Use the original key that worked for storage
            for possible_key in [key, key[:-6] if key.endswith('.ncmat') else None, f"{key}.ncmat"]:
                if possible_key and possible_key in dict.keys(self):
                    super().__setitem__(possible_key, material_dict)
                    break
        return material_dict
    
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
        """Return MaterialSpec-wrapped values."""
        _ensure_initialized()
        for key in super().keys():
            yield self[key]  # Use __getitem__ to get wrapped version
    
    def items(self):
        """Return items with MaterialSpec-wrapped values."""
        _ensure_initialized()
        for key in super().keys():
            yield key, self[key]  # Use __getitem__ to get wrapped version
    
    def get(self, key, default=None):
        _ensure_initialized()
        if key not in self:
            return default
        return self[key]  # Use __getitem__ to get wrapped version
    
    def update(self, *args, **kwargs):
        """Override update to handle initialization properly."""
        _ensure_initialized()
        super().update(*args, **kwargs)

# Add to your module exports or integrate into existing code
__all__ = ['MaterialSpec', 'Hyperlink', 'LazyMaterialsDict']

# Replace materials dict with lazy version
materials = LazyMaterialsDict()