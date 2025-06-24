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
        Process NCMAT files by creating individual templates for each material
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

            # Determine template creation strategy based on section order
            if cell_start < ext_start:
                # @CELL section appears first
                pre_cell_lines = lines[:cell_start + 1]
                post_cell_lines = lines[cell_end:ext_start + 1]
                post_ext_lines = lines[ext_end:] if ext_end else []
                
                self.datatemplate[material] = '\n'.join(
                    pre_cell_lines + 
                    ['**cell_section**'] + 
                    post_cell_lines + 
                    ['**extinction_section**'] + 
                    post_ext_lines
                )
            else:
                # @CUSTOM_CRYSEXTN section appears first
                pre_ext_lines = lines[:ext_start + 1]
                post_ext_lines = lines[ext_end:cell_start + 1]
                post_cell_lines = lines[cell_end:] if cell_end else []
                
                self.datatemplate[material] = '\n'.join(
                    pre_ext_lines + 
                    ['**extinction_section**'] + 
                    post_ext_lines + 
                    ['**cell_section**'] + 
                    post_cell_lines
                )

            # Handle extinction information if present in NCMAT or self.materials
            ext_lines = lines[ext_start + 1:ext_end] if ext_start < len(lines) and ext_start + 1 < ext_end else []
            if ext_lines and ext_lines[0].strip():
                self._extinction_info(material, extinction_lines=ext_lines[0])
            elif any(self.materials[material].get(key) is not None for key in ['ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt']):
                # If extinction parameters are already in self.materials, use them
                self._extinction_info(material)

            # Save original rawdata in nbragg file name
            nc.registerInMemoryFileData(
                self.materials[material]["mat"].replace("ncmat", "nbragg"), 
                self.textdata[material]
            )

            # Apply any user-modified parameters from self.materials
            kwargs = {key: self.materials[material][key] for key in ['a', 'b', 'c', 'ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt'] if self.materials[material][key] is not None}
            if kwargs:
                self._update_ncmat_parameters(material, **kwargs)

    def _update_ncmat_parameters(self, material: str, **kwargs):
        """
        Update the virtual material with lattice and extinction parameters
        
        Args:
            material (str): Name of the material to update
            **kwargs: Additional parameters to update (e.g., a, b, c, ext_l, ext_Gg, ext_L, ext_tilt)
        """
        # Ensure we have a template for this specific material
        if material not in self.datatemplate:
            return

        # Update material parameters if provided in kwargs
        for key in ['ext_l', 'ext_Gg', 'ext_L', 'ext_tilt', 'ext_method', 'a', 'b', 'c']:
            if key in kwargs:
                self.materials[material][key] = float(kwargs[key]) if key != 'ext_method' and key != "ext_tilt" else kwargs[key]

        # Update cell information
        updated_cells = self._cell_info(material, **kwargs)
        
        # Handle extinction information
        updated_ext = self._extinction_info(material, **kwargs) if material in self.extinction or any(kwargs.get(key) for key in ['ext_l', 'ext_Gg', 'ext_L', 'ext_tilt', 'ext_method']) else ""
        
        # Create the updated material text using the material-specific template
        updated_textdata = self.datatemplate[material].replace(
            "**cell_section**", 
            updated_cells
        ).replace(
            "**extinction_section**", 
            updated_ext
        )
        
        # Update the textdata for this specific material
        self.textdata[material] = updated_textdata
        
        # Register the in-memory file with the correct material name
        nc.registerInMemoryFileData(
            self.materials[material]["mat"].replace("ncmat", "nbragg"), 
            updated_textdata
        )

    def _extinction_info(self, material: str, extinction_lines: str = None, **kwargs) -> str:
        """
        Parse and update extinction parameters, storing them directly in self.materials.

        Args:
            material (str): Material name
            extinction_lines (str, optional): Text data from the extinction custom section
            **kwargs: Parameters to update (e.g., ext_l, ext_Gg, ext_L, ext_tilt, ext_method)

        Returns:
            str: Formatted extinction information string
        """
        # Initialize extinction dictionary if not present
        if material not in self.extinction:
            self.extinction[material] = {}

        # Parse extinction lines if provided (from NCMAT file)
        if extinction_lines:
            try:
                method, l, Gg, L, tilt = extinction_lines.strip().split()
                self.extinction[material].update({
                    'method': method,
                    'l': float(l),
                    'Gg': float(Gg),
                    'L': float(L),
                    'tilt': tilt
                })
                # Update self.materials with NCMAT values, preserving any user-provided values
                current_values = {
                    'ext_method': self.materials[material].get('ext_method'),
                    'ext_l': self.materials[material].get('ext_l'),
                    'ext_Gg': self.materials[material].get('ext_Gg'),
                    'ext_L': self.materials[material].get('ext_L'),
                    'ext_tilt': self.materials[material].get('ext_tilt')
                }
                self.materials[material].update({
                    'ext_method': current_values['ext_method'] if current_values['ext_method'] is not None else method,
                    'ext_l': current_values['ext_l'] if current_values['ext_l'] is not None else float(l),
                    'ext_Gg': current_values['ext_Gg'] if current_values['ext_Gg'] is not None else float(Gg),
                    'ext_L': current_values['ext_L'] if current_values['ext_L'] is not None else float(L),
                    'ext_tilt': current_values['ext_tilt'] if current_values['ext_tilt'] is not None else tilt
                })
            except ValueError:
                self.extinction[material] = {}
                # Set defaults in self.materials if no user-provided values exist
                if all(self.materials[material].get(key) is None for key in ['ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt']):
                    self.materials[material].update({
                        'ext_method': None,
                        'ext_l': None,
                        'ext_Gg': None,
                        'ext_L': None,
                        'ext_tilt': None
                    })

        # Update extinction parameters from kwargs (e.g., from user modifications)
        for key, target_key in [('ext_method', 'method'), ('ext_l', 'l'), ('ext_Gg', 'Gg'), ('ext_L', 'L'), ('ext_tilt', 'tilt')]:
            if key in kwargs:
                value = float(kwargs[key]) if key != 'ext_method' and key != 'ext_tilt' else kwargs[key]
                self.extinction[material][target_key] = value
                self.materials[material][key] = value

        # Use self.materials values if available, otherwise apply defaults
        if any(self.materials[material].get(key) is not None for key in ['ext_method', 'ext_l', 'ext_Gg', 'ext_L', 'ext_tilt']):
            self.extinction[material].update({
                'method': self.materials[material].get('ext_method', 'BC_pure'),
                'l': float(self.materials[material].get('ext_l', 0.0)),
                'Gg': float(self.materials[material].get('ext_Gg', 0.0)),
                'L': float(self.materials[material].get('ext_L', 0.0)),
                'tilt': self.materials[material].get('ext_tilt', 'Gauss')
            })
        else:
            # Apply defaults if no values are provided
            defaults = {'method': 'BC_pure', 'l': 0.0, 'Gg': 0.0, 'L': 0.0, 'tilt': 'Gauss'}
            for key, value in defaults.items():
                if key not in self.extinction[material]:
                    self.extinction[material][key] = value
                    self.materials[material][f'ext_{key}'] = value

        # Format the extinction line if extinction parameters are defined
        if self.extinction[material].get('method'):
            method = self.extinction[material]['method']
            l = self.extinction[material]['l']
            Gg = self.extinction[material]['Gg']
            L = self.extinction[material]['L']
            tilt = self.extinction[material]['tilt']
            return f"  {method}  {l:.4f}  {Gg:.4f}  {L:.4f}  {tilt}"
        
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
                     ext_l1, ext_Gg1, ext_L1, ext_tilt1 ... for extinction params
        """
        updated = False
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
            ext_tilt_key = f"ext_tilt{i}"
            
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
                self._update_ncmat_parameters(name,
                                            a=kwargs[lata_key],
                                            b=kwargs[latb_key],
                                            c=kwargs[latc_key])
                updated = True
            elif "a" in kwargs: # for single phase materials
                self._update_ncmat_parameters(name,
                                            a=kwargs["a"],
                                            b=kwargs["b"],
                                            c=kwargs["c"])
                updated = True
            if ext_l_key in kwargs:
                self._update_ncmat_parameters(name,
                                            ext_l=kwargs[ext_l_key],
                                            ext_Gg=kwargs[ext_Gg_key],
                                            ext_L=kwargs[ext_L_key],
                                            ext_tilt=kwargs.get(ext_tilt_key, spec.get('ext_tilt', 'Gauss')))
                updated = True
            elif "ext_l" in kwargs: # for single phase materials
                self._update_ncmat_parameters(name,
                                            ext_l=kwargs["ext_l"],
                                            ext_Gg=kwargs["ext_Gg"],
                                            ext_L=kwargs["ext_L"],
                                            ext_tilt=kwargs.get("ext_tilt", spec.get('ext_tilt', 'Gauss')))
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
        """Plot the cross-section data."""
        import matplotlib.pyplot as plt
        # update lattice and extinction parameters
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

    def _cell_info(self, material: str, **kwargs) -> str:
        """
        Parse crystal cell information and format it for specific output.
        
        Args:
            material (str): Material name
            kwargs (optional): Updates to a, b, c parameters in units of Å
        
        Returns:
            str: Formatted cell information string
        """
        # Load material data directly from NCMAT file to avoid phases_data dependency
        mat_data = nc.load(self.materials[material]["mat"].replace("ncmat", "nbragg"))
        cell_dict = mat_data.info.structure_info
        cell_dict.update(**kwargs)
        return f"  lengths {cell_dict['a']:.4f}  {cell_dict['b']:.4f}  {cell_dict['c']:.4f}  \n  angles {cell_dict['alpha']:.4f}  {cell_dict['beta']:.4f}  {cell_dict['gamma']:.4f}"