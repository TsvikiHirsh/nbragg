from nbragg.utils import materials as materials_dict

import os
import pandas as pd
import numpy as np
import NCrystal as nc
from typing import Dict, Union, Optional, List
from copy import deepcopy

import contextlib
import io
import functools
import sys

def suppress_print(func):
    """Decorator to suppress all stdout/stderr printing from a function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            return func(*args, **kwargs)
    return wrapper

class CrossSection:
    """
    Represents a combination of cross-sections for crystal materials.
    """
    def __init__(self, materials: Union[Dict[str, Union[Dict, dict]], 'CrossSection', None] = None,
                 name: str = None,
                 total_weight: float = 1.,
                 **kwargs):
        """
        Initialize the CrossSection class.
        
        Args:
            materials: Dictionary of material specifications in format:
                {"name": {"mat": material_source, "temp": temp, "mos": mos, "dir1": dir1, "dir2": dir2, "weight": weight}}
                OR {"name": material_dict_from_nbragg_materials}
                OR an instance of the CrossSection class
            name: Name for this cross-section combination.
            **kwargs: Additional materials in format material_name=material_dict_from_nbragg_materials
                      or material_name="material_name_in_nbragg_materials".
        """
        self.name = name
        self.lambda_grid = np.arange(1.0, 10.0, 0.01)  # Default wavelength grid in Ångstroms
        self.mat_data = None  # Single NCrystal scatter object
        self.total_weight = total_weight

        # Initialize materials by combining materials and kwargs
        combined_materials = {}
        
        # Add materials from 'materials' if it is an instance of CrossSection or a dictionary
        if isinstance(materials, CrossSection):
            combined_materials.update(materials.materials)
        elif isinstance(materials, dict):
            combined_materials.update(materials)
        
        # Add materials from kwargs
        for key, value in kwargs.items():
            if isinstance(value, str) and value in materials_dict:
                combined_materials[key] = materials_dict[value]
            else:
                combined_materials[key] = value

        # Process the combined materials dictionary
        self.materials = self._process_materials(combined_materials)
        self.extinction = {}

        # Create virtual material
        self._create_virtual_materials()
        
        # Initialize weights
        self.weights = pd.Series(dtype=float)
        self._set_weights()
        self._generate_cfg_string()
        self._load_material_data()
        self._populate_material_data()

    def _process_materials(self, materials: Dict[str, Union[Dict, dict]]) -> Dict[str, Dict]:
        """Process materials dictionary while preserving relative weights."""
        processed = {}
        raw_total_weight = 0
        
        # First pass: process specifications without normalizing weights
        for name, spec in materials.items():
            if isinstance(spec, dict) and not spec.get('mat'):
                processed[name] = {
                    'mat': spec.get('mat'),
                    'temp': spec.get('temp', 300.),
                    'mos': spec.get('mos', None),
                    'dir1': spec.get('dir1', None),
                    'dir2': spec.get('dir2', None),
                    'dirtol': spec.get('dirtol', None),
                    'theta': spec.get('theta', None),
                    'phi': spec.get('phi', None),
                    'a': spec.get('a', None),
                    'b': spec.get('b', None),
                    'c': spec.get('c', None),
                    'ext_method': spec.get('ext_method', None),
                    'ext_l': spec.get('ext_l', None),
                    'ext_Gg': spec.get('ext_Gg', None),
                    'ext_L': spec.get('ext_L', None),
                    'ext_tilt': spec.get('ext_tilt', None),
                    'weight': spec.get('weight', 1.0)
                }
                raw_total_weight += processed[name]['weight']
            elif isinstance(spec, CrossSection):
                for material_name, material_spec in spec.materials.items():
                    processed[f"{name}_{material_name}"] = material_spec.copy()
                    raw_total_weight += material_spec['weight']
            else:
                if not isinstance(spec, dict):
                    raise ValueError(f"Material specification for {name} must be a dictionary")
                
                material = spec.get('mat')
                if isinstance(material, dict):
                    material = material.get('mat')
                elif isinstance(material, str):
                    material = self._resolve_material(material)
                    
                weight = float(spec.get('weight', 1.0))
                raw_total_weight += weight
                
                processed[name] = {
                    'mat': material,
                    'temp': spec.get('temp', 300.),
                    'mos': spec.get('mos', None),
                    'dir1': spec.get('dir1', None),
                    'dir2': spec.get('dir2', None),
                    'dirtol': spec.get('dirtol', None),
                    'theta': spec.get('theta', None),
                    'phi': spec.get('phi', None),
                    'a': spec.get('a', None),
                    'b': spec.get('b', None),
                    'c': spec.get('c', None),
                    'ext_method': spec.get('ext_method', None),
                    'ext_l': spec.get('ext_l', None),
                    'ext_Gg': spec.get('ext_Gg', None),
                    'ext_L': spec.get('ext_L', None),
                    'ext_tilt': spec.get('ext_tilt', None),
                    'weight': weight
                }
        
        # Second pass: normalize weights while preserving relative proportions
        if raw_total_weight > 0:
            for spec in processed.values():
                spec['weight'] = (spec['weight'] / raw_total_weight)

        return processed
    
    def _create_virtual_materials(self):
        """
        Process NCMAT files by creating individual templates for each material, preserving or updating @CUSTOM_CRYSEXTN sections.
        """
        # Initialize dictionaries to store material-specific data
        self.textdata = {}
        self.datatemplate = {}
        
        for material in self.materials:
            # Save entire input text
            self.textdata[material] = nc.createTextData(self.materials[material]["mat"]).rawData
            
            # Split input into lines
            lines = self.textdata[material].split('\n')
            
            # Find @CELL section
            cell_start = None
            cell_end = None
            for i, line in enumerate(lines):
                if line.strip().startswith('@CELL'):
                    cell_start = i
                elif cell_start is not None and line.strip().startswith('@'):
                    cell_end = i
                    break
            
            # Find @CUSTOM_CRYSEXTN section
            ext_start = None
            ext_end = None
            for i, line in enumerate(lines):
                if line.strip().startswith('@CUSTOM_CRYSEXTN'):
                    ext_start = i
                elif ext_start is not None and line.strip().startswith('@'):
                    ext_end = i
                    break
            
            # Handle cases where sections might be missing
            ext_start = len(lines) if ext_start is None else ext_start
            ext_end = len(lines) if ext_end is None else ext_end
            cell_start = len(lines) if cell_start is None else cell_start
            cell_end = len(lines) if cell_end is None else cell_end

            # Create template based on section order
            if cell_start < ext_start:
                # @CELL section appears before @CUSTOM_CRYSEXTN or no extinction section
                pre_cell_lines = lines[:cell_start + 1]
                post_cell_lines = lines[cell_end:ext_start]
                ext_lines = lines[ext_start:ext_end] if ext_start < len(lines) else []
                post_ext_lines = lines[ext_end:] if ext_end < len(lines) else []
                
                self.datatemplate[material] = '\n'.join(
                    pre_cell_lines + 
                    ['**cell_section**'] + 
                    post_cell_lines + 
                    (['@CUSTOM_CRYSEXTN', '**extinction_section**'] if ext_lines else ['**extinction_section**']) + 
                    post_ext_lines
                )
            else:
                # @CUSTOM_CRYSEXTN section appears before @CELL or no cell section
                pre_ext_lines = lines[:ext_start]
                ext_lines = lines[ext_start:ext_end] if ext_start < len(lines) else []
                post_ext_lines = lines[ext_end:cell_start]
                post_cell_lines = lines[cell_end:] if cell_end < len(lines) else []
                
                self.datatemplate[material] = '\n'.join(
                    pre_ext_lines + 
                    (['@CUSTOM_CRYSEXTN', '**extinction_section**'] if ext_lines else ['**extinction_section**']) + 
                    post_ext_lines + 
                    ['**cell_section**'] + 
                    post_cell_lines
                )

            # Handle extinction information
            ext_lines_content = lines[ext_start + 1:ext_end] if ext_start < len(lines) and ext_start + 1 < ext_end and lines[ext_start + 1].strip() else []
            has_ext_params = any(self.materials[material].get(key) is not None for key in ['ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt'])
            
            if ext_lines_content:
                # Existing @CUSTOM_CRYSEXTN section in NCMAT
                self._extinction_info(material, extinction_lines=ext_lines_content[0])
            elif has_ext_params:
                # No @CUSTOM_CRYSEXTN in NCMAT, but ext_ parameters provided
                self._extinction_info(material)
            else:
                # No extinction information; apply defaults if needed later
                self.extinction[material] = {}

            # Save original rawdata in nbragg file name
            nc.registerInMemoryFileData(
                self.materials[material]["mat"].replace("ncmat", "nbragg"), 
                self.textdata[material]
            )

            # Apply any user-modified parameters from self.materials
            kwargs = {key: self.materials[material][key] for key in ['a', 'b', 'c', 'ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt'] if self.materials[material][key] is not None}
            if kwargs:
                self._update_ncmat_parameters(material, **kwargs)

    def update(self):
        """
        Update the CrossSection object after modifying self.materials.
        Reprocesses material parameters, updates virtual materials, and reloads data.
        """
        # Update virtual materials with current parameters
        for material in self.materials:
            kwargs = {key: self.materials[material][key] for key in ['a', 'b', 'c', 'ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt'] if self.materials[material].get(key) is not None}
            if kwargs:
                self._update_ncmat_parameters(material, **kwargs)
        
        # Update weights and data
        self._set_weights()
        self._generate_cfg_string()
        self._load_material_data()
        self._populate_material_data()

    @suppress_print
    def _update_ncmat_parameters(self, material: str, **kwargs):
        """
        Update the virtual material with lattice and extinction parameters.

        Args:
            material (str): Name of the material to update
            **kwargs: Additional parameters to update (e.g., a, b, c, ext_l, ext_Gg, ext_L, ext_tilt)
        """
        # Ensure we have a template for this specific material
        if material not in self.datatemplate:
            return

        # Update material parameters if provided in kwargs
        for key in ['ext_l', 'ext_Gg', 'ext_L', 'ext_tilt', 'ext_method', 'a', 'b', 'c']:
            if key in kwargs and kwargs[key] is not None:
                self.materials[material][key] = float(kwargs[key]) if key not in ['ext_method', 'ext_tilt'] else kwargs[key]

        # Update cell information
        updated_cells = self._cell_info(material, **kwargs)
        
        # Handle extinction information
        updated_ext = self._extinction_info(material, **kwargs) if any(kwargs.get(key) for key in ['ext_l', 'ext_Gg', 'ext_L', 'ext_tilt', 'ext_method']) or material in self.extinction else ""
        
        # Create the updated material text using the material-specific template
        updated_textdata = self.datatemplate[material].replace(
            "**cell_section**", 
            updated_cells
        ).replace(
            "@CUSTOM_CRYSEXTN\n**extinction_section**" if "@CUSTOM_CRYSEXTN" in self.datatemplate[material] else "**extinction_section**", 
            "@CUSTOM_CRYSEXTN\n" + updated_ext if updated_ext else ""
        )
        
        # Update the textdata for this specific material
        self.textdata[material] = updated_textdata
        
        # Register the in-memory file with the correct material name
        nc.registerInMemoryFileData(
            self.materials[material]["mat"].replace("ncmat", "nbragg"), 
            updated_textdata
        )

    @suppress_print
    def _extinction_info(self, material: str, extinction_lines: str = None, **kwargs) -> str:
        """
        Parse and update extinction parameters, storing them directly in self.materials.
        Extinction is only valid if the material has crystallographic structure info.

        @CUSTOM_CRYSEXTN formats (from CrysExtn plugin source code):
        - Sabine_uncorr  l  G  L  rect/tri
        - Sabine_corr    l  g  L
        - BC_mix         l  g  L  Gauss/Lorentz/Fresnel
        - BC_pure        l  g  L  Gauss/Lorentz/Fresnel
        - BC_mod         l  g  L  Gauss/Lorentz/Fresnel
        - CR             l  g  L
        - RED_orig       l  R  T
        - RED            l  R  T  c

        Where:
        - l: crystallite size/path length (Å)
        - g/G/R: mosaicity parameter (1/rad) or equivalent
        - L/T: grain size/path length (Å)
        - c: correlation parameter for RED (0 to 1)
        - rect/tri: tilt distribution (for Sabine_uncorr only)
        - Gauss/Lorentz/Fresnel: orientation method (for BC models only)
        """
        mat_path = self.materials[material]["mat"]
        mat_data = nc.load(mat_path)

        # Check if material has structure info
        try:
            _ = mat_data.info.structure_info
            has_structure = True
        except Exception:
            has_structure = False

        # Initialize extinction dictionary if not present
        if material not in self.extinction:
            self.extinction[material] = {}

        # Define default extinction parameters
        defaults = {
            'ext_method': 'BC_pure',
            'ext_l': 2500.0,      # crystallite size/path length (Å)
            'ext_g': 150.0,       # mosaicity parameter (1/rad)
            'ext_L': 100000.0,    # grain size/path length (Å)
            'ext_dist': 'Gauss',  # distribution: rect/tri for Sabine_uncorr, Gauss/Lorentz/Fresnel for BC
            'ext_c': 0.5          # correlation for RED
        }

        # Supported extinction methods and distributions
        supported_methods = ['Sabine_uncorr', 'Sabine_corr', 'BC_mix', 'BC_pure', 'BC_mod', 'CR', 'RED_orig', 'RED']
        methods_with_dist = ['Sabine_uncorr', 'BC_mix', 'BC_pure', 'BC_mod']
        methods_with_four_params = ['Sabine_corr', 'CR', 'RED_orig']
        sabine_distributions = ['rect', 'tri']
        bc_distributions = ['Gauss', 'Lorentz', 'Fresnel']

        # Parse extinction lines from NCMAT (if present)
        if extinction_lines:
            try:
                parts = extinction_lines.strip().split()
                if len(parts) < 4:
                    raise ValueError(f"Expected at least 4 parameters, got {len(parts)}")
                
                method = parts[0]
                if method not in supported_methods:
                    method = defaults['ext_method']
                
                l = float(parts[1]) if float(parts[1]) > 0 else defaults['ext_l']
                g = float(parts[2]) if float(parts[2]) > 0 else defaults['ext_g']
                L = float(parts[3]) if float(parts[3]) > 0 else defaults['ext_L']
                
                dist = None
                c = None
                if method in methods_with_dist:
                    if len(parts) >= 5:
                        dist = parts[4]
                        if method == 'Sabine_uncorr' and dist not in sabine_distributions:
                            dist = 'rect'
                        elif method in ['BC_mix', 'BC_pure', 'BC_mod'] and dist not in bc_distributions:
                            dist = 'Gauss'
                    else:
                        dist = 'rect' if method == 'Sabine_uncorr' else 'Gauss'
                elif method == 'RED' and len(parts) >= 5:
                    c = float(parts[4]) if 0 <= float(parts[4]) <= 1 else defaults['ext_c']
                
                self.extinction[material].update({
                    'method': method,
                    'l': l,
                    'g': g,
                    'L': L,
                    'dist': dist if method in methods_with_dist else None,
                    'c': c if method == 'RED' else None
                })
            except (ValueError, IndexError) as e:
                import warnings
                warnings.warn(f"Could not parse extinction line for {material}: {e}")
                self.extinction[material] = {}

        # First pass: process ext_method if provided to ensure it's available for validation
        if 'ext_method' in kwargs and kwargs['ext_method'] is not None:
            method_value = kwargs['ext_method'] if kwargs['ext_method'] in supported_methods else defaults['ext_method']
            self.extinction[material]['method'] = method_value
            self.materials[material]['ext_method'] = method_value
        
        # Apply user-provided kwargs, prioritizing over NCMAT values
        for key, target_key in [
            ('ext_l', 'l'),
            ('ext_g', 'g'),
            ('ext_Gg', 'g'),  # Support both ext_g and ext_Gg
            ('ext_L', 'L'),
            ('ext_dist', 'dist'),
            ('ext_tilt', 'dist'),  # Support ext_tilt for backward compatibility
            ('ext_c', 'c')
        ]:
            if key in kwargs and kwargs[key] is not None:
                value = kwargs[key]
                if key in ['ext_dist', 'ext_tilt']:
                    # Get current method from various sources (including what we just set above)
                    current_method = (
                        self.extinction[material].get('method') or 
                        self.materials[material].get('ext_method') or 
                        defaults['ext_method']
                    )
                    if current_method == 'Sabine_uncorr':
                        value = value if value in sabine_distributions else 'rect'
                    elif current_method in ['BC_mix', 'BC_pure', 'BC_mod']:
                        value = value if value in bc_distributions else 'Gauss'
                    else:
                        value = None
                elif key == 'ext_c':
                    try:
                        value = float(value)
                        value = value if 0 <= value <= 1 else defaults['ext_c']
                    except (ValueError, TypeError):
                        value = defaults['ext_c']
                elif key in ['ext_l', 'ext_g', 'ext_Gg', 'ext_L']:
                    try:
                        value = float(value)
                        value = value if value > 0 else defaults[key if key != 'ext_Gg' else 'ext_g']
                    except (ValueError, TypeError):
                        value = defaults[key if key != 'ext_Gg' else 'ext_g']
                
                self.extinction[material][target_key] = value
                # Store in materials with the key name used (ext_g for ext_Gg, ext_dist for ext_tilt)
                if key == 'ext_Gg':
                    self.materials[material]['ext_g'] = value
                elif key == 'ext_tilt':
                    self.materials[material]['ext_dist'] = value
                    self.materials[material]['ext_tilt'] = value  # Keep both for backward compatibility
                else:
                    self.materials[material][key] = value

        # Apply defaults if extinction parameters are partially specified
        if has_structure and any(self.materials[material].get(k) is not None for k in defaults.keys()):
            for k, d in defaults.items():
                if self.materials[material].get(k) is None:
                    self.materials[material][k] = d
            
            method = self.materials[material]['ext_method']
            if method in methods_with_dist and self.materials[material].get('ext_dist') is None:
                default_dist = 'rect' if method == 'Sabine_uncorr' else 'Gauss'
                self.materials[material]['ext_dist'] = default_dist
                self.materials[material]['ext_tilt'] = default_dist  # Keep both
            if method == 'RED' and self.materials[material].get('ext_c') is None:
                self.materials[material]['ext_c'] = defaults['ext_c']
            
            self.extinction[material].update({
                'method': method,
                'l': float(self.materials[material]['ext_l']),
                'g': float(self.materials[material].get('ext_g', self.materials[material].get('ext_Gg', defaults['ext_g']))),
                'L': float(self.materials[material]['ext_L']),
                'dist': self.materials[material].get('ext_dist') if method in methods_with_dist else None,
                'c': self.materials[material].get('ext_c') if method == 'RED' else None
            })

        # Return formatted extinction block only if the material is crystalline
        if not has_structure:
            return ""

        if self.extinction[material].get('method'):
            method = self.extinction[material]['method']
            line = (
                f"  {method}  "
                f"{self.extinction[material]['l']:.4f}  "
                f"{self.extinction[material]['g']:.4f}  "
                f"{self.extinction[material]['L']:.4f}"
            )
            if method in methods_with_dist and self.extinction[material].get('dist'):
                line += f"  {self.extinction[material]['dist']}"
            elif method == 'RED' and self.extinction[material].get('c') is not None:
                line += f"  {self.extinction[material]['c']:.4f}"
            return line
        return ""

    def _resolve_material(self, material: str) -> str:
        """Resolve material specification to filename."""
        if material.endswith('.ncmat'):
            return material
            
        mat_info = self._get_material_info(material)
        if mat_info:
            return mat_info.get('mat')
        return material

    def _set_weights(self):
        """Set weights from processed materials."""
        if not self.materials:
            self.weights = pd.Series(dtype=float)
            return
            
        # Apply total_weight to all material weights
        self.weights = pd.Series({name: spec['weight'] * self.total_weight
                                for name, spec in self.materials.items()})
        
        # Update the material weights
        for name, weight in self.weights.items():
            self.materials[name]['weight'] = weight

    def __add__(self, other: 'CrossSection') -> 'CrossSection':
        """Add two CrossSection objects."""
        combined_materials = {}
        
        # Add materials from both objects
        for name, spec in self.materials.items():
            combined_materials[name] = deepcopy(spec)
            
        # Add materials from other, ensuring unique names
        for name, spec in other.materials.items():
            new_name = name
            counter = 1
            while new_name in combined_materials:
                new_name = f"{name}_{counter}"
                counter += 1
            combined_materials[new_name] = deepcopy(spec)
        
        return CrossSection(combined_materials)

    def __mul__(self, scalar: float) -> 'CrossSection':
        """Multiply CrossSection by a scalar."""
        new_materials = deepcopy(self.materials)
        result = CrossSection(new_materials, total_weight=scalar)
        return result
    
    def __rmul__(self, scalar) -> 'CrossSection':
        # For commutative multiplication (scalar * material)
        return self.__mul__(scalar)
    
    def _generate_cfg_string(self):
        """
        Generate configuration strings using NCrystal phase notation with consistent phase ordering.
        Stores individual phase configurations in self.phases dictionary and
        creates a combined configuration string in self.cfg_string.
        """
        if not self.materials:
            self.cfg_string = ""
            self.phases = {}
            return

        # Sort materials by their keys to ensure consistent ordering
        sorted_materials = dict(sorted(self.materials.items()))
        
        phase_parts = []
        self.phases = {}
        # Calculate the sum of weights for normalization
        total = sum(spec['weight'] for spec in sorted_materials.values())
        
        for name, spec in sorted_materials.items():
            material = spec['mat']
            if not material:
                continue
                
            # Normalize the weight for NCrystal configuration
            normalized_weight = spec['weight'] / total if total > 0 else spec['weight']
            phase = f"{normalized_weight}*{material}"
            single_phase = f"{material}"
            
            # Collect material-specific parameters
            params = []
            if spec['temp'] is not None:
                params.append(f"temp={spec['temp']}K")
                
            # Determine if the material is oriented
            mos = spec.get('mos', None)
            dir1 = spec.get('dir1', None)
            dir2 = spec.get('dir2', None)
            dirtol = spec.get('dirtol', None)
            theta = spec.get('theta', None)
            phi = spec.get('phi', None)
            
            is_oriented = mos is not None or dir1 is not None or dir2 is not None
            if is_oriented:
                # Apply default values if not provided
                mos = mos if mos is not None else 0.001
                dir1 = dir1 if dir1 is not None else (0, 0, 1)
                dir2 = dir2 if dir2 is not None else (1, 0, 0)
                dirtol = dirtol if dirtol is not None else 1.
                theta = theta if theta is not None else 0.
                phi = phi if phi is not None else 0.
                
                # Format the orientation vectors with NCrystal-specific notation
                orientation = self.format_orientations(dir1, dir2, theta=theta, phi=phi)
                dir1_str = f"@crys_hkl:{orientation['dir1'][0]:.8f},{orientation['dir1'][1]:.8f},{orientation['dir1'][2]:.8f}@lab:0,0,1"
                dir2_str = f"@crys_hkl:{orientation['dir2'][0]:.8f},{orientation['dir2'][1]:.8f},{orientation['dir2'][2]:.8f}@lab:0,1,0"
                
                params.append(f"mos={mos}deg")
                params.append(f"dirtol={dirtol}deg")
                params.append(f"dir1={dir1_str}")
                params.append(f"dir2={dir2_str}")
                
            # Combine parameters with the phase if any exist
            if params:
                phase += f";{';'.join(sorted(params))}"  # Sort parameters for consistency
                single_phase += f";{';'.join(sorted(params))}"
                
            # Store the individual phase configuration in the dictionary and replace materials with virtual mat
            self.phases[name] = single_phase.replace("ncmat", "nbragg")
            # Add to the list for the combined configuration string
            phase_parts.append(phase)
            
        # Generate the complete configuration string
        self.cfg_string = f"phases<{'&'.join(phase_parts)}>" if phase_parts else ""
        # replace materials with virtual materials
        self.cfg_string = self.cfg_string.replace("ncmat", "nbragg")

    @suppress_print
    def _load_material_data(self):
        """Load the material data using NCrystal with the phase configuration."""
        if self.cfg_string:
            self.mat_data = nc.load(self.cfg_string)

    @suppress_print
    def _populate_material_data(self):
        """Populate cross section data using NCrystal phases."""
        if not self.cfg_string:
            self.table = pd.DataFrame(index=self.lambda_grid)
            self.table.index.name = "wavelength"
            return

        xs = {}

        # Load all phases in the final weights order
        self.phases_data = {name: nc.load(self.phases[name]) for name in self.weights.index}
        for phase in self.weights.index:
            xs[phase] = self._calculate_cross_section(self.lambda_grid, self.phases_data[phase])

        # Calculate total
        xs["total"] = self._calculate_cross_section(self.lambda_grid, self.mat_data)

        # Build DataFrame in the correct order from the start
        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"

        if not hasattr(self, "atomic_density"):
            self.atomic_density = self.mat_data.info.factor_macroscopic_xs

    def _calculate_cross_section(self, wl, mat):
        """Calculate cross-section using NCrystal's xsect method."""
        xs = mat.scatter.xsect(wl=wl, direction=(0,0,1)) + mat.absorption.xsect(wl=wl, direction=(0,0,1))
        return np.nan_to_num(xs,0.)

    @suppress_print
    def __call__(self, wl: np.ndarray, **kwargs):
        """
        Update configuration if parameters change and return cross-section.
        
        Args:
            wl: Wavelength array
            **kwargs: Material-specific parameters in format:
                    η1, η2, ... for mosaic spread of materials 1, 2, ...
                    θ1, θ2, ... for theta values of materials 1, 2, ...
                    ϕ1, ϕ2, ... for phi values of materials 1, 2, ...
                    temp1, temp2, ... for temperatures of materials 1, 2, ...
                    a1, a2, ... for lattice parameter of materials 1, 2 ...
                    ext_l1, ext_Gg1, ext_L1, ext_tilt1, ext_method1 ... for extinction params
        """
        updated = False
        # Check for parameter updates
        material_names = list(self.materials.keys())
        for i, name in enumerate(material_names, 1):
            spec = self.materials[name]
            
            # Check for material-specific parameters
            temp_key = f"temp"  # all phase temperatures are updated to the same value
            mos_key = f"η{i}"
            theta_key = f"θ{i}"
            phi_key = f"ϕ{i}"
            lata_key = f"a{i}"
            latb_key = f"b{i}"
            latc_key = f"c{i}"
            ext_l_key = f"ext_l{i}"
            ext_Gg_key = f"ext_Gg{i}"
            ext_L_key = f"ext_L{i}"
            ext_tilt_key = f"ext_tilt{i}"
            ext_method_key = f"ext_method{i}"
            
            # Update temperature
            if temp_key in kwargs and kwargs[temp_key] != spec['temp']:
                spec['temp'] = kwargs[temp_key]
                updated = True
            
            # Update mosaic spread
            if mos_key in kwargs and kwargs[mos_key] != spec['mos']:
                spec['mos'] = kwargs[mos_key]
                updated = True
            
            # Update theta
            if theta_key in kwargs and kwargs[theta_key] != spec['theta']:
                spec['theta'] = kwargs[theta_key]
                updated = True
            
            # Update phi
            if phi_key in kwargs and kwargs[phi_key] != spec['phi']:
                spec['phi'] = kwargs[phi_key]
                updated = True
            
            # Update phase weight
            phase_name = name.replace("-", "")
            if phase_name in kwargs and kwargs[phase_name] != spec["weight"]:
                spec['weight'] = kwargs[phase_name]
                updated = True
            
            # Update lattice parameters
            if lata_key in kwargs:
                self._update_ncmat_parameters(name,
                                            a=kwargs[lata_key],
                                            b=kwargs[latb_key],
                                            c=kwargs[latc_key])
                updated = True
            elif "a" in kwargs:  # for single phase materials
                self._update_ncmat_parameters(name,
                                            a=kwargs["a"],
                                            b=kwargs["b"],
                                            c=kwargs["c"])
                updated = True
            
            # Update extinction parameters
            if ext_l_key in kwargs or ext_Gg_key in kwargs or ext_L_key in kwargs or ext_tilt_key in kwargs or ext_method_key in kwargs:
                ext_params = {}
                if ext_l_key in kwargs:
                    ext_params['ext_l'] = kwargs[ext_l_key]
                elif ext_l_key not in kwargs and ext_Gg_key in kwargs:
                    ext_params['ext_l'] = spec.get('ext_l', 2500.0)  # Default from _extinction_info
                if ext_Gg_key in kwargs:
                    ext_params['ext_Gg'] = kwargs[ext_Gg_key]
                if ext_L_key in kwargs:
                    ext_params['ext_L'] = kwargs[ext_L_key]
                elif ext_L_key not in kwargs and (ext_l_key in kwargs or ext_Gg_key in kwargs):
                    ext_params['ext_L'] = spec.get('ext_L', 100000.0)  # Default from _extinction_info
                if ext_tilt_key in kwargs:
                    ext_params['ext_tilt'] = kwargs[ext_tilt_key]
                if ext_method_key in kwargs:
                    ext_params['ext_method'] = kwargs[ext_method_key]
                self._update_ncmat_parameters(name, **ext_params)
                updated = True
            elif "ext_l" in kwargs:  # for single phase materials
                ext_params = {}
                if "ext_l" in kwargs:
                    ext_params['ext_l'] = kwargs["ext_l"]
                if "ext_Gg" in kwargs:
                    ext_params['ext_Gg'] = kwargs["ext_Gg"]
                if "ext_L" in kwargs:
                    ext_params['ext_L'] = kwargs["ext_L"]
                else:
                    ext_params['ext_L'] = spec.get('ext_L', 100000.0)  # Default
                if "ext_tilt" in kwargs:
                    ext_params['ext_tilt'] = kwargs["ext_tilt"]
                if "ext_method" in kwargs:
                    ext_params['ext_method'] = kwargs["ext_method"]
                self._update_ncmat_parameters(name, **ext_params)
                updated = True

        if updated:
            self._set_weights()
            self._generate_cfg_string()
            self._load_material_data()
            self._populate_material_data()

        return self._calculate_cross_section(wl, self.mat_data)
    
    @staticmethod
    def _get_material_info(material_key: str) -> Dict:
        """Get material information from the materials dictionary."""
        material_info = materials_dict.get(material_key)
        
        if not material_info:
            for info in materials_dict.values():
                if (info.get('formula') == material_key or 
                    info.get('name') == material_key or 
                    info.get('mat') == material_key):
                    return info
        
        return material_info

    def plot(self, **kwargs):
        """
        Plot the weighted neutron cross-section data for each phase and the total.

        This method will:
        1. Update lattice and extinction parameters for each material (if applicable).
        2. Load and populate the cross-section data table.
        3. Plot each phase's weighted cross-section in the same order as the table columns.
        4. Plot the total cross-section as a thicker dark line.
        5. Generate a legend with phase names and their weight percentages.

        Parameters
        ----------
        title : str, optional
            Title of the plot. Defaults to the cross-section object's `name`.
        ylabel : str, optional
            Y-axis label. Defaults to ``"σ [barn]"``.
        xlabel : str, optional
            X-axis label. Defaults to ``"Wavelength [Å]"``.
        lw : float, optional
            Base line width for phase curves. Defaults to ``1.0``.
        **kwargs
            Additional keyword arguments passed to ``pandas.DataFrame.plot`` for the phase curves.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib Axes object containing the plot.

        Notes
        -----
        - The order of phases in the plot and legend is preserved according to the
        order of the columns in ``self.table`` (excluding the "total" column).
        - The "total" curve is always plotted last, with a distinct style and label.

        Examples
        --------
        >>> xs = 0.0275 * nbragg.CrossSection(celloluse=nbragg.materials["Cellulose_C6O5H10.ncmat"]) \
        ...     + (1 - 0.00275) * nbragg.CrossSection(α=nbragg.materials["Fe_sg229_Iron-alpha_CrysExtn1.ncmat"])
        >>> ax = xs.plot(title="Cross-section comparison", lw=1.5)
        >>> ax.figure.show()
        """
        import matplotlib.pyplot as plt

        # Update parameters if possible
        try:
            for material in self.materials:
                self._update_ncmat_parameters(material)
        except Exception:
            pass

        self._load_material_data()
        self._populate_material_data()

        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "σ [barn]")
        xlabel = kwargs.pop("xlabel", "Wavelength [Å]")
        lw = kwargs.pop("lw", 1.0)

        fig, ax = plt.subplots()

        # Ensure weights are aligned with table column order (excluding total)
        phase_cols = [col for col in self.table.columns if col != "total"]
        weights_aligned = self.weights.reindex(phase_cols)

        # Plot each phase component
        if phase_cols:
            self.table[phase_cols].mul(weights_aligned).plot(ax=ax, lw=lw, **kwargs)

        # Plot total
        self.table["total"].plot(ax=ax, color="0.2", lw=lw * 1.2, label="Total")

        # Axis labels and title
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        # Legend
        legend_labels = [f"{mat}: {weights_aligned[mat] * 100:.3f}%" for mat in phase_cols] + ["Total"]
        ax.legend(legend_labels)

        return ax


    @staticmethod
    def _normalize_vector(vector: Union[List[float], List[str]]) -> List[float]:
        """Normalizes a vector to have a length of 1.
        
        Args:
            vector: List of numbers (as floats or strings)
            
        Returns:
            List[float]: Normalized vector
        """
        # Convert strings to floats if necessary
        vec_float = [float(x) if isinstance(x, str) else x for x in vector]
        magnitude = sum(x**2 for x in vec_float) ** 0.5
        if magnitude == 0:
            return [0.0, 0.0, 0.0]
        return [x / magnitude for x in vec_float]

    def _rotate_vector(self, vec: List[float], phi: float = 0.0, theta: float = 0.0) -> List[float]:
        """Rotates a vector by angles phi (around z-axis) and theta (around y-axis)."""
        # Ensure vector components are floats
        vec = [float(x) if isinstance(x, str) else x for x in vec]
        
        # Convert angles from degrees to radians
        phi = np.radians(float(phi))
        theta = np.radians(float(theta))
        
        # Rotation matrix around z-axis
        Rz = np.array([
            [np.cos(phi), -np.sin(phi), 0],
            [np.sin(phi),  np.cos(phi), 0],
            [0,            0,           1]
        ])
        
        # Rotation matrix around y-axis
        Ry = np.array([
            [ np.cos(theta), 0, np.sin(theta)],
            [ 0,             1, 0            ],
            [-np.sin(theta), 0, np.cos(theta)]
        ])
        
        # Apply rotations: first around z, then around y
        rotated_vec = Ry @ (Rz @ np.array(vec, dtype=float))
        return rotated_vec.tolist()

    def format_orientations(self, dir1: Union[List[float], List[str]] = None, 
                            dir2: Union[List[float], List[str]] = None,
                            phi: Union[float, str] = 0.0, 
                            theta: Union[float, str] = 0.0) -> Dict[str, List[float]]:
        """Converts dir1 and dir2 vectors to NCrystal orientation format with optional rotation."""
        if dir1 is None:
            dir1 = [0.0, 0.0, 1.0]
        if dir2 is None:
            dir2 = [1.0, 0.0, 0.0]

        # Convert any string values to floats and normalize
        dir1 = self._normalize_vector([float(x) if isinstance(x, str) else x for x in dir1])
        dir2 = self._normalize_vector([float(x) if isinstance(x, str) else x for x in dir2])
        phi = float(phi) if isinstance(phi, str) else phi
        theta = float(theta) if isinstance(theta, str) else theta

        # Apply rotations if specified
        if phi != 0 or theta != 0:
            dir1 = self._rotate_vector(dir1, phi, theta)
            dir2 = self._rotate_vector(dir2, phi, theta)

        # Return vectors without any string formatting for easy processing
        return {
            'dir1': dir1,
            'dir2': dir2
        }

    @classmethod
    def _normalize_mtex_vector(cls, vector):
        """Normalize a vector to unit length."""
        vec = np.array(vector)
        magnitude = np.linalg.norm(vec)
        return (vec / magnitude).tolist() if magnitude > 0 else vec.tolist()

    @classmethod
    def from_mtex(cls, csv_file, material, powder_phase=True, short_name=None):
        """
        Create a CrossSection from MTEX CSV orientation data.
        
        Parameters:
        -----------
        csv_file : str
            Path to the CSV file containing orientation components
        material : dict
            Base material dictionary with existing properties
        powder_phase : bool, optional
            Whether to add a non-oriented powder phase with complementary weight (default True)
        short_name : str, optional
            Short name for the phase (e.g., 'γ' for gamma)
        
        Returns:
        --------
        CrossSection
            CrossSection object with materials from CSV data
        """
        # Read the CSV file
        try:
            df = pd.read_csv(csv_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"CSV file not found: {csv_file}")
        
        # Handle column name variations
        column_mapping = {
            'alpha_mtex': ['alpha_mtex', 'alpha'],
            'beta_mtex': ['beta_mtex', 'beta'],
            'gamma_mtex': ['gamma_mtex', 'gamma'],
            'volume_mtex': ['volume_mtex', 'volume'],
            'fwhm': ['fwhm', 'fwhm_mtex'],
            'xh': ['xh'],
            'xk': ['xk'],
            'xl': ['xl'],
            'yh': ['yh'],
            'yk': ['yk'],
            'yl': ['yl']
        }
        
        # Find the correct column names
        def find_column(key_list):
            for key in key_list:
                if key in df.columns:
                    return key
            raise KeyError(f"Could not find column for {key_list}")
        
        # Map columns
        try:
            alpha_col = find_column(column_mapping['alpha_mtex'])
            beta_col = find_column(column_mapping['beta_mtex'])
            gamma_col = find_column(column_mapping['gamma_mtex'])
            volume_col = find_column(column_mapping['volume_mtex'])
            fwhm_col = find_column(column_mapping['fwhm'])
            xh_col = find_column(column_mapping['xh'])
            xk_col = find_column(column_mapping['xk'])
            xl_col = find_column(column_mapping['xl'])
            yh_col = find_column(column_mapping['yh'])
            yk_col = find_column(column_mapping['yk'])
            yl_col = find_column(column_mapping['yl'])
        except KeyError:
            # If specific orientation columns are not found, return a CrossSection with base material
            return cls({short_name or material['name']: material}, name=short_name)
        
        # Normalize volumes to ensure they sum to 1 or less
        total_volume = df[volume_col].sum()
        if total_volume > 1:
            df[volume_col] = df[volume_col] / total_volume
        
        # Prepare materials dictionary
        materials = {}
        
        # Process each row for oriented phases
        for i, row in df.iterrows():
            # Create a copy of the base material
            updated_material = material.copy()
            
            # Extract weight
            weight = row[volume_col]
            
            # MTEX to NCrystal coordinate transformation
            dir1 = cls._normalize_mtex_vector([row[xh_col], row[xk_col], row[xl_col]])
            dir2 = cls._normalize_mtex_vector([row[yh_col], row[yk_col], row[yl_col]])
            
            # Create material name
            material_name = f"{short_name or material['name']}{i}"
            
            # Update material dictionary with parameters
            updated_material.update({
                'mat': material.get('mat', ''),  # Use existing mat or empty string
                'temp': material.get('temp', 300),  # Default to 300 if not specified
                'mos': row[fwhm_col],
                'dir1': dir1,
                'dir2': dir2,
                'alpha': row[alpha_col],
                'beta': row[beta_col],
                'gamma': row[gamma_col],
                'dirtol': 1.0,
                'theta': 0.0,
                'phi': 0.0,
                'weight': weight
            })
            
            # Add to materials dictionary
            materials[material_name] = updated_material
        
        # Add non-oriented powder phase if requested
        if powder_phase:
            background_weight = 1.0 - total_volume if total_volume <= 1 else 0.0
            if background_weight > 0:
                background_material = material.copy()
                background_material.update({
                    'mat': material.get('mat', ''),
                    'temp': material.get('temp', 300),
                    'mos': None,  # No mosaicity for powder phase
                    'dir1': None,  # No orientation
                    'dir2': None,
                    'alpha': None,
                    'beta': None,
                    'gamma': None,
                    'dirtol': None,
                    'theta': None,
                    'phi': None,
                    'weight': background_weight
                })
                materials[f"{short_name or material['name']}_powder"] = background_material
        
        # Return CrossSection with materials
        return cls(materials, name=short_name)

    @classmethod
    def _estimate_mosaicity(cls, df):
        """Estimate mosaicity from the dataframe."""
        # If FWHM column exists, use it directly
        fwhm_cols = ['fwhm', 'fwhm_mtex']
        for col in fwhm_cols:
            if col in df.columns:
                return df[col].mean()
        
        # If no FWHM, try to estimate from volume spread
        if 'volume' in df.columns:
            volume_std = df['volume'].std()
            base_mosaicity = 5.0  # degrees, adjust as needed
            adjusted_mosaicity = base_mosaicity * (1 + volume_std * 10)
            return min(adjusted_mosaicity, 50.0)
        
        # If no information available, return None
        return None
    
    def _cell_info(self, material, **kwargs):
        """
        Retrieve crystallographic cell information if available,
        otherwise return an empty string.
        """
        mat_path = self.materials[material]["mat"]

        # Always load original NCMAT file (not the virtual nbragg)
        mat_data = nc.load(mat_path)

        # Try to get structure info
        try:
            cell_dict = mat_data.info.structure_info
        except Exception:
            # No structure info available (e.g., amorphous or molecular material)
            return ""

        # Apply overrides from kwargs
        cell_dict.update(**kwargs)

        return (
            f"  lengths {cell_dict['a']:.4f}  {cell_dict['b']:.4f}  {cell_dict['c']:.4f}  \n"
            f"  angles {cell_dict['alpha']:.4f}  {cell_dict['beta']:.4f}  {cell_dict['gamma']:.4f}"
        )
