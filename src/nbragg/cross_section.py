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
                # Replace the string with the actual material dictionary
                combined_materials[key] = materials_dict[value]
            else:
                # Assume the provided value is already a valid material dictionary
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
        """Process and normalize the materials dictionary."""
        processed = {}
        total_weight = 0
        
        # First pass: process specifications and sum weights
        for name, spec in materials.items():
            if isinstance(spec, dict) and not spec.get('mat'):
                # Handle direct nbragg.materials dictionary input
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
                total_weight += processed[name]['weight']
            elif isinstance(spec, CrossSection):
                # Handle CrossSection instance
                for material_name, material_spec in spec.materials.items():
                    processed[f"{name}_{material_name}"] = material_spec
                total_weight += sum(material_spec['weight'] for material_spec in spec.materials.values())
            else:
                if not isinstance(spec, dict):
                    raise ValueError(f"Material specification for {name} must be a dictionary")
                
                material = spec.get('mat')
                if isinstance(material, dict):
                    # Handle nbragg.materials object in 'mat' key
                    material = material.get('mat')
                elif isinstance(material, str):
                    material = self._resolve_material(material)
                    
                weight = float(spec.get('weight', 1.0))
                total_weight += weight
                
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
        
        # Second pass: normalize weights
        if total_weight > 0:
            for spec in processed.values():
                spec['weight'] /= total_weight
                
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
            
        self.weights = pd.Series({name: spec['weight'] 
                                for name, spec in self.materials.items()})

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
        
        # Create new instance with combined materials
        return CrossSection(combined_materials, name=f"{self.name}+{other.name}")

    def __mul__(self, scalar: float) -> 'CrossSection':
        """Multiply CrossSection by a scalar."""
        new_materials = {}
        
        for name, spec in self.materials.items():
            new_spec = deepcopy(spec)
            new_spec['weight'] *= scalar
            new_materials[name] = new_spec
        
        return CrossSection(new_materials, name=f"{scalar}*{self.name}")

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

        for name, spec in self.materials.items():
            material = spec['mat']
            if not material:
                continue

            # Build the base phase configuration
            phase = f"{spec['weight']}*{material}"
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
                


                # Format the orientations
                orientation = self.format_orientations(dir1, dir2,theta=theta,phi=phi)
                dir1_str = orientation['dir1']
                dir2_str = orientation['dir2']
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
        return mat.scatter.xsect(wl=wl, direction=(0,0,1)) + mat.absorption.xsect(wl=wl, direction=(0,0,1))

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
            if theta_key in kwargs and kwargs[theta] != spec['theta']:
                spec['theta'] = kwargs[theta_key]
                updated = True
            if phi_key in kwargs and kwargs[phi] != spec['phi']:
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


    def _normalize_vector(self,vec):
        """Normalizes a vector."""
        norm = np.linalg.norm(vec)
        if norm == 0:
            return vec
        return [x / norm for x in vec]

    def from_maud(self,maud_line=None):
        """Parses a MAUD line or returns default vectors aligned with the beam.
        
        Args:
            maud_line (str, optional): MAUD string. Defaults to a crystal aligned with the beam.
            
        Returns:
            dict: A dictionary with normalized dir1 and dir2 vectors.
        """
        default_maud = (
            "Vol:0.10, EA_ZXZ:(0.00 0.00 0.00), "
            "x||(1.0000 0.0000 0.0000), y||(0.0000 1.0000 0.0000), z||(0.0000 0.0000 1.0000)"
        )
        
        if maud_line is None:
            maud_line = default_maud
        
        try:
            # Extract the vectors from the MAUD line
            parts = maud_line.split(',')
            x_vector = [float(x) for x in parts[2].split('||')[1].strip('()').split()]
            z_vector = [float(x) for x in parts[4].split('||')[1].strip('()').split()]
            
            # Normalize the vectors
            dir1 = self._normalize_vector(z_vector)
            dir2 = self._normalize_vector(x_vector)
            
            return {'dir1': dir1, 'dir2': dir2}
        
        except (IndexError, ValueError):
            raise ValueError(
                "Invalid MAUD line format. Expected format example: "
                "'Vol:0.10, EA_ZXZ:(77.21 45.31 268.14), x||(1.9669 0.7061 2.0107), y||(-0.5429 2.8119 -0.4564), z||(-2.0607 -0.0669 2.0394)'"
            )

    def _rotate_vector(self,vec, phi=0., theta=0.):
        """Rotates a vector by angles phi (around z-axis) and theta (around y-axis)."""
        # Convert angles from degrees to radians
        phi = np.radians(phi)
        theta = np.radians(theta)
        
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
        rotated_vec = Ry @ (Rz @ np.array(vec))
        return rotated_vec.tolist()

    def format_orientations(self, dir1=None, dir2=None, phi=0, theta=0):
        """Converts dir1 and dir2 vectors to NCrystal orientation format with optional rotation.
        
        Args:
            dir1 (list of float, optional): Normalized vector for dir1. Defaults to z-axis aligned.
            dir2 (list of float, optional): Normalized vector for dir2. Defaults to x-axis perpendicular.
            phi (float, optional): Rotation around the z-axis in degrees. Defaults to 0.
            theta (float, optional): Rotation around the y-axis in degrees. Defaults to 0.
            
        Returns:
            dict: A dictionary with the NCrystal formatted dir1 and dir2 strings.
        """
        if dir1 is None:
            dir1 = [0.0, 0.0, 1.0]
        if dir2 is None:
            dir2 = [1.0, 0.0, 0.0]
        
        # Rotate the vectors
        dir1_rotated = self._rotate_vector(dir1, phi, theta)
        dir2_rotated = self._rotate_vector(dir2, phi, theta)
        
        # Format the vectors
        dir1_str = ",".join(f"{coord:.4f}" for coord in dir1_rotated)
        dir2_str = ",".join(f"{coord:.4f}" for coord in dir2_rotated)
        
        return {
            'dir1': f"@crys_hkl:{dir1_str}@lab:0,0,1",
            'dir2': f"@crys_hkl:{dir2_str}@lab:0,1,0"
        }