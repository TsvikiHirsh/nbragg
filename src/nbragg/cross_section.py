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
                    'weight': weight
                }
        
        # Second pass: normalize weights while preserving relative proportions
        if raw_total_weight > 0:
            for spec in processed.values():
                spec['weight'] = (spec['weight'] / raw_total_weight)
                
        return processed

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
    
    def _generate_cfg_string(self):
        """
        Generate configuration strings using NCrystal phase notation.
        Stores individual phase configurations in self.phases dictionary and
        creates a combined configuration string in self.cfg_string.
        """
        if not self.materials:
            self.cfg_string = ""
            self.phases = {}
            return

        phase_parts = []
        self.phases = {}

        # Calculate the sum of weights for normalization
        total = sum(spec['weight'] for spec in self.materials.values())

        for name, spec in self.materials.items():
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
                params.append(f"temp={int(spec['temp'])}K")

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
                dir1_str = f"@crys_hkl:{orientation['dir1'][0]:.4f},{orientation['dir1'][1]:.4f},{orientation['dir1'][2]:.4f}@lab:0,0,1"
                dir2_str = f"@crys_hkl:{orientation['dir2'][0]:.4f},{orientation['dir2'][1]:.4f},{orientation['dir2'][2]:.4f}@lab:0,1,0"

                params.append(f"mos={mos}deg")
                params.append(f"dirtol={dirtol}deg")
                params.append(f"dir1={dir1_str}")
                params.append(f"dir2={dir2_str}")

            # Combine parameters with the phase if any exist
            if params:
                phase += f";{';'.join(params)}"
                single_phase += f";{';'.join(params)}"

            # Store the individual phase configuration in the dictionary
            self.phases[name] = single_phase

            # Add to the list for the combined configuration string
            phase_parts.append(phase)

        # Generate the complete configuration string
        self.cfg_string = f"phases<{'&'.join(phase_parts)}>" if phase_parts else ""


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
        
        mat = nc.load(self.cfg_string)
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
        
        self.atomic_density = mat.info.factor_macroscopic_xs

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
        
        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "$\sigma$ [barn]")
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

    def _transform_lab_coordinates(self, vector: List[float]) -> List[float]:
        """Transform from HIPPO lab coordinates to NCrystal lab coordinates."""
        # Ensure vector components are floats
        vector = [float(x) if isinstance(x, str) else x for x in vector]
        
        transform = np.array([
            [0, 0, 1],  # HIPPO's z becomes NCrystal's x
            [0, 1, 0],  # HIPPO's y stays as NCrystal's y
            [1, 0, 0]   # HIPPO's x becomes NCrystal's z
        ])
        
        return (transform @ np.array(vector, dtype=float)).tolist()

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

        # Transform to NCrystal lab coordinates
        dir1_lab = self._transform_lab_coordinates(dir1)
        dir2_lab = self._transform_lab_coordinates(dir2)

        # Return vectors without any string formatting for easy processing
        return {
            'dir1': dir1_lab,
            'dir2': dir2_lab
        }


    def from_maud(self, maud_line: str = None, dirtol: float = None, 
                 mos: float = None, suffix: str = None) -> 'CrossSection':
        """Parses a MAUD line and converts orientations from HIPPO to NCrystal convention."""
        default_maud = (
            "Vol:0.10, EA_ZXZ:(0.00 0.00 0.00), "
            "x||(1.0000 0.0000 0.0000), y||(0.0000 1.0000 0.0000), z||(0.0000 0.0000 1.0000)"
        )
        
        if maud_line is None:
            maud_line = default_maud
        
        try:
            # Extract volume fraction
            parts = maud_line.split(',')
            vol_str = parts[0].split(':')[1].strip()
            vol = float(vol_str)
            
            # Extract vectors from MAUD line and convert to floats
            x_vector = [float(x) for x in parts[2].split('||')[1].strip('()').split()]
            y_vector = [float(x) for x in parts[3].split('||')[1].strip('()').split()]
            
            # Get orientation strings
            orientations = self.format_orientations(
                dir1=x_vector,  # HIPPO x (beam) direction
                dir2=y_vector   # HIPPO y (vertical) direction
            )
            
            # Update materials
            updated_materials = {}
            for name, spec in self.materials.items():
                updated_spec = spec.copy()
                updated_spec['dir1'] = orientations['dir1']
                updated_spec['dir2'] = orientations['dir2']
                if dirtol is not None:
                    updated_spec['dirtol'] = float(dirtol)
                if mos is not None:
                    updated_spec['mos'] = float(mos)
                
                new_name = f"{name}{suffix}" if suffix else name
                updated_materials[new_name] = updated_spec
            
            return CrossSection(updated_materials, total_weight=vol)
            
        except (IndexError, ValueError) as e:
            raise ValueError(
                f"Invalid MAUD line format: {str(e)}. Expected format example: "
                "'Vol:0.10, EA_ZXZ:(77.21 45.31 268.14), x||(1.9669 0.7061 2.0107), "
                "y||(-0.5429 2.8119 -0.4564), z||(-2.0607 -0.0669 2.0394)'"
            )
        
    def from_mtex(self, csv_file, material, short_name=None):
        """
        Extract orientation and additional information from a MTEX CSV file.
        
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
        dict
            Updated material dictionary with additional orientation information
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
            # If specific orientation columns are not found, return the material as-is
            return {short_name or material['name']: material}
        
        # Normalize volumes to ensure they sum to 1 or less
        total_volume = df[volume_col].sum()
        if total_volume > 1:
            print(f"Warning: Total volume exceeds 1.0. Normalizing...")
            df[volume_col] = df[volume_col] / total_volume
        
        # Use the last row for orientation information
        last_row = df.iloc[-1]
        
        # Extract Euler angles and volume
        euler_angles = [last_row[alpha_col], last_row[beta_col], last_row[gamma_col]]
        weight = last_row[volume_col]
        
        # Prepare the material dictionary
        updated_material = material.copy()
        
        # Add additional keys
        updated_material.update({
            'mos': self._estimate_mosaicity(df),
            'dir1': self._extract_vector(last_row, 'x'),
            'dir2': self._extract_vector(last_row, 'y'),
            'dirtol': None,  # You may want to add logic to extract this if available
            'theta': euler_angles[1],  # Beta angle as theta
            'phi': euler_angles[2],    # Gamma angle as phi
            'weight': weight
        })
        
        # Use short name or material name as the key
        return {short_name or material['name']: updated_material}

    def _extract_vector(self, row, axis):
        """
        Extract vector components for a given axis.
        
        Parameters:
        -----------
        row : pandas.Series
            DataFrame row containing vector components
        axis : str
            Axis to extract ('x', 'y', or 'z')
        
        Returns:
        --------
        list
            Vector components [h, k, l]
        """
        vector_cols = [
            f'{axis}h_mtex', f'{axis}h', 
            f'{axis}k_mtex', f'{axis}k', 
            f'{axis}l_mtex', f'{axis}l'
        ]
        
        # Try to find vector components
        for col in vector_cols:
            if col in row.index:
                # If the column exists and has a non-None value
                if pd.notna(row[col]):
                    return [
                        row.get(f'{axis}h_mtex', row.get(f'{axis}h', 1.0 if axis == 'x' else 0.0)),
                        row.get(f'{axis}k_mtex', row.get(f'{axis}k', 0.0 if axis == 'y' else 1.0)),
                        row.get(f'{axis}l_mtex', row.get(f'{axis}l', 0.0 if axis == 'z' else 1.0))
                    ]
        
        # Default vector if no specific components found
        default_vectors = {
            'x': [1.0, 0.0, 0.0],
            'y': [0.0, 1.0, 0.0],
            'z': [0.0, 0.0, 1.0]
        }
        return default_vectors[axis]

    def _estimate_mosaicity(self, df):
        """
        Estimate mosaicity from the dataframe.
        
        Parameters:
        -----------
        df : pandas.DataFrame
            Dataframe containing orientation components
        
        Returns:
        --------
        float or None
            Estimated mosaicity value
        """
        # If FWHM column exists, use it directly
        fwhm_cols = ['fwhm', 'fwhm_mtex']
        for col in fwhm_cols:
            if col in df.columns:
                return df[col].mean()
        
        # If no FWHM, return None
        return None