from nbragg.utils import materials as materials_dict

import os
import pandas as pd
import numpy as np
import NCrystal as nc
from typing import Dict, Union, Optional, List
from copy import deepcopy

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

        # create virtual material
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
                    'a': spec.get('a',None),
                    'b': spec.get('b',None),
                    'c': spec.get('c',None),
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
                    'a': spec.get('a',None),
                    'b': spec.get('b',None),
                    'c': spec.get('c',None),
                    'weight': weight
                }


        
        # Second pass: normalize weights while preserving relative proportions
        if raw_total_weight > 0:
            for spec in processed.values():
                spec['weight'] = (spec['weight'] / raw_total_weight)

        return processed
    

    def _create_virtual_materials(self):
        """
        Process NCMAT files by creating individual templates for each material.
        Handles both crystalline (with @CELL and/or @CUSTOM_CRYSEXTN) and 
        non-crystalline (no @CELL or @CUSTOM_CRYSEXTN) materials.
        """
        # Initialize dictionaries to store material-specific data
        self.textdata = {}
        self.datatemplate = {}
        
        for material in self.materials:
            # Save entire input text
            self.textdata[material] = nc.createTextData(self.materials[material]["mat"]).rawData
            
            # Split input into lines
            lines = self.textdata[material].split('\n')
            
            # Check for presence of @CELL and @CUSTOM_CRYSEXTN sections
            has_cell = any(line.strip().startswith('@CELL') for line in lines)
            has_extinction = any(line.strip().startswith('@CUSTOM_CRYSEXTN') for line in lines)
            
            if has_cell or has_extinction:
                # Handle crystalline materials with @CELL and/or @CUSTOM_CRYSEXTN
                cell_start = None
                cell_end = None
                ext_start = None
                ext_end = None
                
                # Find section boundaries
                for i, line in enumerate(lines):
                    if line.strip().startswith('@CELL'):
                        cell_start = i
                    elif cell_start is not None and line.strip().startswith('@') and i > cell_start:
                        cell_end = i
                        break
                    if line.strip().startswith('@CUSTOM_CRYSEXTN'):
                        ext_start = i
                    elif ext_start is not None and line.strip().startswith('@') and i > ext_start:
                        ext_end = i
                        break
                
                # Default to end of file if sections are not terminated
                cell_start = cell_start if cell_start is not None else len(lines)
                cell_end = cell_end if cell_end is not None else len(lines)
                ext_start = ext_start if ext_start is not None else len(lines)
                ext_end = ext_end if ext_end is not None else len(lines)

                # Determine template creation strategy based on section order
                if has_cell and (not has_extinction or cell_start < ext_start):
                    # @CELL section appears first or is the only section
                    pre_cell_lines = lines[:cell_start + 1]
                    post_cell_lines = lines[cell_end:ext_start + 1] if has_extinction else lines[cell_end:]
                    post_ext_lines = lines[ext_end:] if has_extinction else []
                    
                    self.datatemplate[material] = '\n'.join(
                        pre_cell_lines + 
                        ['**cell_section**'] + 
                        post_cell_lines + 
                        (['**extinction_section**'] + post_ext_lines if has_extinction else [])
                    )
                elif has_extinction:
                    # @CUSTOM_CRYSEXTN section appears first or is the only section
                    pre_ext_lines = lines[:ext_start + 1]
                    post_ext_lines = lines[ext_end:cell_start + 1] if has_cell else lines[ext_end:]
                    post_cell_lines = lines[cell_end:] if has_cell else []
                    
                    self.datatemplate[material] = '\n'.join(
                        pre_ext_lines + 
                        ['**extinction_section**'] + 
                        post_ext_lines + 
                        (['**cell_section**'] + post_cell_lines if has_cell else [])
                    )

                # Handle extinction information if present
                if has_extinction:
                    ext_lines = lines[ext_start + 1] if ext_start + 1 < len(lines) else ""
                    if ext_lines:
                        self._extinction_info(material, extinction_lines=ext_lines)
            else:
                # Non-crystalline material: no @CELL or @CUSTOM_CRYSEXTN
                # Use the raw text as the template, no sections to replace
                self.datatemplate[material] = self.textdata[material]
                # No extinction data to process

            # Register the in-memory file with the virtual material name
            nc.registerInMemoryFileData(
                self.materials[material]["mat"].replace("ncmat", "nbragg"), 
                self.textdata[material]
            )

    def _update_ncmat_parameters(self, material: str, **kwargs):
        """
        Update the virtual material with parameters. For non-crystalline materials,
        skip lattice and extinction updates as they are not applicable.
        
        Args:
            material (str): Name of the material to update
            **kwargs: Additional parameters (e.g., a, b, c for lattice, l, Gg, L for extinction)
        """
        # Ensure we have a template for this specific material
        if material not in self.datatemplate:
            return

        # Check if the material is crystalline (has placeholders for sections)
        is_crystalline = '**cell_section**' in self.datatemplate[material] or '**extinction_section**' in self.datatemplate[material]
        
        if is_crystalline:
            # Update cell information if lattice parameters are provided
            if 'a' in kwargs or 'b' in kwargs or 'c' in kwargs:
                updated_cells = self._cell_info(material, **kwargs)
            else:
                # Use existing cell section if no update is provided
                cell_start = self.textdata[material].find('**cell_section**')
                if cell_start != -1:
                    updated_cells = '**cell_section**'
                else:
                    updated_cells = ""

            # Handle extinction information if present and updated
            if material in self.extinction and any(k in kwargs for k in ['l', 'Gg', 'L']):
                updated_ext = self._extinction_info(material, **kwargs)
            else:
                # Use existing extinction section if no update is provided
                ext_start = self.textdata[material].find('**extinction_section**')
                if ext_start != -1:
                    updated_ext = '**extinction_section**'
                else:
                    updated_ext = ""
            
            # Create the updated material text using the material-specific template
            updated_textdata = self.datatemplate[material].replace(
                "**cell_section**", 
                updated_cells
            ).replace(
                "**extinction_section**", 
                updated_ext
            )
        else:
            # Non-crystalline material: no updates to lattice or extinction, use original text
            updated_textdata = self.datatemplate[material]

        # Update the textdata for this specific material
        self.textdata[material] = updated_textdata
        
        # Register the in-memory file with the correct material name
        nc.registerInMemoryFileData(
            self.materials[material]["mat"].replace("ncmat", "nbragg"), 
            updated_textdata
        )

    def _extinction_info(self, material: str, extinction_lines: str = None, **kwargs) -> str:
        """
        Parse and update the extinction lines. Returns empty string for non-crystalline materials
        unless explicitly updated.

        Args:
            material (str): Material name
            extinction_lines (str, optional): Text data from the extinction custom section
            **kwargs: Additional parameters (e.g., l, Gg, L)
        """
        # Initialize extinction data if not already present
        if material not in self.extinction:
            self.extinction[material] = {}

        if extinction_lines:
            # Parse extinction data if provided (crystalline material)
            method, l, Gg, L, tilt = extinction_lines.split()
            self.extinction[material] = dict(method=method, l=float(l), Gg=float(Gg), L=float(L), tilt=tilt)
        elif not self.extinction[material] and not kwargs:
            # No extinction data and no updates: non-crystalline material
            return ""
        
        # Update extinction parameters if provided
        if kwargs and any(k in kwargs for k in ['l', 'Gg', 'L']):
            self.extinction[material].update(**{k: v for k, v in kwargs.items() if k in ['l', 'Gg', 'L']})
            # Assume default values for missing parameters if extinction is being added
            self.extinction[material].setdefault('method', 'gaussian')
            self.extinction[material].setdefault('l', 0.0)
            self.extinction[material].setdefault('Gg', 0.0)
            self.extinction[material].setdefault('L', 0.0)
            self.extinction[material].setdefault('tilt', '0')

        # Return formatted extinction string if data exists
        if self.extinction[material]:
            method = self.extinction[material]["method"]
            l = self.extinction[material]["l"]
            Gg = self.extinction[material]["Gg"]
            L = self.extinction[material]["L"]
            tilt = self.extinction[material]["tilt"]
            return f"  {method}  {l}  {Gg}  {L}  {tilt}"
        return ""

    def _cell_info(self, material: str, **kwargs) -> str:
        """
        Generate or update the @CELL section string. Only applicable to crystalline materials.
        """
        # Check if the material is crystalline
        if '**cell_section**' not in self.datatemplate[material]:
            return ""  # Non-crystalline materials don’t have a cell section
        
        # Default lattice parameters from kwargs or existing material spec
        a = kwargs.get('a', self.materials[material].get('a', 1.0))
        b = kwargs.get('b', self.materials[material].get('b', 1.0))
        c = kwargs.get('c', self.materials[material].get('c', 1.0))
        
        # For simplicity, assume cubic lattice unless specified otherwise in the future
        return f"  lengths {a} {b} {c}\n  angles 90 90 90"
        

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

    def _load_material_data(self):
        """Load the material data using NCrystal with the phase configuration."""
        if self.cfg_string:
            self.mat_data = nc.load(self.cfg_string)

    def _populate_material_data(self):
        """Populate cross section data using NCrystal phases."""
        if not self.cfg_string:
            self.table = pd.DataFrame(index=self.lambda_grid)
            self.table.index.name = "wavelength"
            return

        xs = {}

        # Process each phase separately
        self.phases_data = {name:nc.load(self.phases[name]) for name in self.phases} 
        for phase in self.phases:
            xs[phase] = self._calculate_cross_section(self.lambda_grid, self.phases_data[phase])
        
        # Calculate total
        xs["total"] = self._calculate_cross_section(self.lambda_grid, self.mat_data)
        
        # Create DataFrame with all phases
        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"
        
        if len(self.table.columns) > 1:
            self.table.columns = self.weights.index.to_list() + ["total"]
        else:
            self.table.columns = ["total"]
        
        if not hasattr(self,"atomic_density"):
            self.atomic_density = self.mat_data.info.factor_macroscopic_xs

    def _calculate_cross_section(self, wl, mat):
        """Calculate cross-section using NCrystal's xsect method."""
        xs = mat.scatter.xsect(wl=wl, direction=(0,0,1)) + mat.absorption.xsect(wl=wl, direction=(0,0,1))
        return np.nan_to_num(xs,0.)

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
                     ext_l1, ext_Gg1, ext_L1 ... for extinction params
        """
        updated = False
        direction = None
        # Check for parameter updates
        material_names = list(self.materials.keys())
        for i, name in enumerate(material_names, 1):
            spec = self.materials[name]
            
            # Check for material-specific parameters
            temp_key = f"temp" # all phase temperatures are updated to the same value
            mos_key = f"η{i}"
            theta_key = f"θ{i}"
            phi_key = f"ϕ{i}"
            lata_key = f"a{i}"
            latb_key = f"b{i}"
            latc_key = f"c{i}"
            ext_l_key = f"ext_l{i}"
            ext_Gg_key = f"ext_Gg{i}"
            ext_L_key = f"ext_L{i}"
            
            if temp_key in kwargs and kwargs[temp_key] != spec['temp']:
                spec['temp'] = kwargs[temp_key]
                updated = True
            if mos_key in kwargs and kwargs[mos_key] != spec['mos']:
                spec['mos'] = kwargs[mos_key]
                updated = True
            if theta_key in kwargs and kwargs[theta_key] != spec['theta']:
                spec['theta'] = kwargs[theta_key]
                updated = True
            if phi_key in kwargs and kwargs[phi_key] != spec['phi']:
                spec['phi'] = kwargs[phi_key]
                updated = True
            phase_name = name.replace("-", "")
            if phase_name in kwargs and kwargs[phase_name] != spec["weight"]:
                spec['weight'] = kwargs[phase_name]
                updated = True
            if lata_key in kwargs:
                self._update_ncmat_parameters(name,a=kwargs[lata_key],b=kwargs[latb_key],c=kwargs[latc_key])
                updated = True
            elif "a" in kwargs: # for single phase materials
                self._update_ncmat_parameters(name,a=kwargs["a"],b=kwargs["b"],c=kwargs["c"])
                updated = True
            if ext_l_key in kwargs:
                self._update_ncmat_parameters(name,l=kwargs[ext_l_key],Gg=kwargs[ext_Gg_key],L=kwargs[ext_L_key])
                updated = True
            elif "ext_l" in kwargs: # for single phase materials
                self._update_ncmat_parameters(name,l=kwargs["ext_l"],Gg=kwargs["ext_Gg"],L=kwargs["ext_L"])
                updated = True

        if updated:
            self._generate_cfg_string()
            self._load_material_data()
            

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
        """Plot the cross-section data."""
        import matplotlib.pyplot as plt
        # update lattice parameters
        try:
            for material in self.materials:
                self._update_ncmat_parameters(material)
        except:
            pass
        self._load_material_data()
        self._populate_material_data()
        
        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "σ [barn]")
        xlabel = kwargs.pop("xlabel", "Wavelength [Å]")
        lw = kwargs.pop("lw", 1.)

        fig, ax = plt.subplots()

        # Plot each material component with reduced line width
        if len(self.table.columns) > 1:
            self.table.iloc[:, :-1].mul(self.weights).plot(ax=ax, lw=lw, **kwargs)

        # Plot the total curve with a thicker line width and distinct color
        self.table["total"].plot(ax=ax, color="0.2", lw=lw*1.2, label="Total")

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        legend_labels = [f"{material}: {weight*100:.1f}%" 
                        for material, weight in self.weights.items()] + ["Total"]
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
    def from_mtex(cls, csv_file, material, short_name=None):
        """
        Create a CrossSection from MTEX CSV orientation data.
        
        Parameters:
        -----------
        csv_file : str
            Path to the CSV file containing orientation components
        material : dict
            Base material dictionary with existing properties
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
            'volume_mtex': ['volume_mtex', 'volume']
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
        except KeyError:
            # If specific orientation columns are not found, return a CrossSection with base material
            return cls({short_name or material['name']: material}, name=short_name)
        
        # Normalize volumes to ensure they sum to 1 or less
        total_volume = df[volume_col].sum()
        if total_volume > 1:
            df[volume_col] = df[volume_col] / total_volume
        
        # Prepare materials dictionary
        materials = {}
        
        # Estimate overall mosaicity for reference
        overall_mosaicity = cls._estimate_mosaicity(df)
        
        # Process each row
        for i, row in df.iterrows():
            # Create a copy of the base material
            updated_material = material.copy()
            
            # Extract weight
            weight = row[volume_col]
            
            # Estimate mosaicity for this specific row (or use overall if not possible)
            mos = cls._estimate_mosaicity(df.loc[[i]]) or overall_mosaicity
            
            # MTEX to NCrystal coordinate transformation
            # h normal becomes beam direction (z)
            # k normal becomes x direction 
            # l normal becomes y direction
            dir1 = cls._normalize_mtex_vector([row.get('zh', 0), row.get('zk', 0), row.get('zl', 1)])
            dir2 = cls._normalize_mtex_vector([row.get('yh', 0), row.get('yk', 1), row.get('yl', 0)])
            
            # Create material name
            material_name = f"{short_name or material['name']}{i+1}"
            
            # Update material dictionary
            updated_material.update({
                'mos': mos,
                'dir1': dir1,
                'dir2': dir2,
                'dirtol': None,
                'theta': None,
                'phi': None,
                'weight': weight
            })
            
            # Add to materials dictionary
            materials[material_name] = updated_material
        
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