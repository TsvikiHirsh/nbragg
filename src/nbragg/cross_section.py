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
    def __init__(self, materials: Union[Dict[str, Union[Dict, dict]], None] = None,
                name: str = "",
                temp: float = 300.0,
                mos: Union[float, Dict[str, float], None] = None,
                k: Union[float, Dict[str, float], None] = None,
                l: Union[float, Dict[str, float], None] = None,
                dirtol: float = 1.,
                **kwargs):
        """
        Initialize the CrossSection class.
        
        Args:
            materials: Dictionary of material specifications in format:
                {"name": {"mat": material_source, "temp": temp, "mos": mos, "k": k, "l": l, "weight": weight}}
                OR {"name": material_dict_from_nbragg_materials}
            name: Name for this cross section combination
            temp: Default temperature if not specified per material
            mos, k, l: Default crystal orientation parameters if not specified per material
            dirtol: Direction tolerance in degrees
            **kwargs: Additional materials in format material_name=material_dict_from_nbragg_materials
        """
        self.name = name
        self.default_temp = temp
        self.default_mos = mos
        self.default_k = k
        self.default_l = l
        self.dirtol = dirtol
        self.L = 9.  # TODO: replace this hack
        
        self.lambda_grid = np.arange(1.0, 10.0, 0.01)  # Default wavelength grid in Ångstroms
        self.matdata = None  # Single NCrystal scatter object
        
        # Combine materials from dict and kwargs
        combined_materials = {}
        if materials:
            if isinstance(materials, dict):
                combined_materials.update(materials)
        combined_materials.update(kwargs)
        
        # Process materials dictionary
        self.materials = self._process_materials(combined_materials or {})
        
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
                    'temp': spec.get('temp', self.default_temp),
                    'mos': spec.get('mos', self.default_mos),
                    'k': spec.get('k', self.default_k),
                    'l': spec.get('l', self.default_l),
                    'weight': spec.get('weight', 1.0)
                }
                total_weight += processed[name]['weight']
            else:
                if not isinstance(spec, dict):
                    raise ValueError(f"Material specification for {name} must be a dictionary")
                
                material = spec.get('mat')
                if isinstance(material, dict):
                    # Handle nbragg.materials object in 'mat' key
                    material = material.get('filename')
                elif isinstance(material, str):
                    material = self._resolve_material(material)
                    
                weight = float(spec.get('weight', 1.0))
                total_weight += weight
                
                processed[name] = {
                    'mat': material,
                    'temp': spec.get('temp', self.default_temp),
                    'mos': spec.get('mos', self.default_mos),
                    'k': spec.get('k', self.default_k),
                    'l': spec.get('l', self.default_l),
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
            return mat_info.get('filename')
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
        """Generate a single configuration string using NCrystal phase notation."""
        if not self.materials:
            self.cfg_string = ""
            return
            
        phase_parts = []
        
        for name, spec in self.materials.items():
            material = spec['mat']
            if not material:
                continue
                
            phase = f"{spec['weight']}*{material}"
            
            # Add material-specific parameters
            params = []
            
            if spec['temp'] is not None:
                params.append(f"temp={spec['temp']}K")
            
            if all(param is not None for param in (spec['mos'], spec['k'], spec['l'])):
                params.append(f"mos={spec['mos']}deg")
                params.append(f"dirtol={self.dirtol}deg")
                params.append(f"dir1=@crys_hkl:0,{spec['k']},{spec['l']}@lab:0,0,1")
                params.append(f"dir2=@crys_hkl:0,-1,1@lab:0,1,0")
            
            if params:
                phase += f";{';'.join(params)}"
            
            phase_parts.append(phase)
        
        self.cfg_string = f"phases<{'&'.join(phase_parts)}>" if phase_parts else ""
        self.phases = phase_parts

    def _load_material_data(self):
        """Load the material data using NCrystal with the phase configuration."""
        if self.cfg_string:
            self.matdata = nc.load(self.cfg_string)
            self.phases = {name:{"mat":nc.NCMATComposer.from_info(phase[1]).load(),"weight":phase[0]} 
                         for name,phase in zip(self.phases,self.matdata.info.phases)}

    def _populate_material_data(self):
        """Populate cross section data using NCrystal phases."""
        if not self.cfg_string:
            self.table = pd.DataFrame(index=self.lambda_grid)
            self.table.index.name = "wavelength"
            return

        mat = nc.load(self.cfg_string)
        xs = {}

        # Process each phase separately
        for phase, phase_info in self.phases.items():
            xs[phase] = self._calculate_cross_section(self.lambda_grid, phase_info["mat"])

        # Calculate total
        xs["total"] = self._calculate_cross_section(self.lambda_grid, self.matdata)

        # Create DataFrame with all phases
        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"
        if len(self.table.columns) > 1:
            self.table.columns = self.weights.index.to_list() + ["total"]
        else:
            self.table.columns = ["total"]

        self.atomic_density = mat.info.factor_macroscopic_xs

    def _calculate_cross_section(self, wl, mat, direction=None):
        """Calculate cross-section using NCrystal's xsect method."""
        return mat.scatter.xsect(wl=wl, direction=direction) + mat.absorption.xsect(wl=wl, direction=direction)

    def __call__(self, wl: np.ndarray, **kwargs):
        """
        Update configuration if parameters change and return cross-section.
        
        Args:
            wl: Wavelength array
            **kwargs: Material-specific parameters in format:
                     mos1, mos2, ... for mosaic spread of materials 1, 2, ...
                     k1, k2, ... for k values of materials 1, 2, ...
                     l1, l2, ... for l values of materials 1, 2, ...
                     temp1, temp2, ... for temperatures of materials 1, 2, ...
        """
        updated = False
        direction = None
        
        # Check for parameter updates
        material_names = list(self.materials.keys())
        for i, name in enumerate(material_names, 1):
            spec = self.materials[name]
            
            # Check for material-specific parameters
            temp_key = f"temp{i}"
            mos_key = f"mos{i}"
            k_key = f"k{i}"
            l_key = f"l{i}"
            
            if temp_key in kwargs and kwargs[temp_key] != spec['temp']:
                spec['temp'] = kwargs[temp_key]
                updated = True
            if mos_key in kwargs and kwargs[mos_key] != spec['mos']:
                spec['mos'] = kwargs[mos_key]
                updated = True
            if k_key in kwargs and kwargs[k_key] != spec['k']:
                spec['k'] = kwargs[k_key]
                updated = True
            if l_key in kwargs and kwargs[l_key] != spec['l']:
                spec['l'] = kwargs[l_key]
                updated = True

        if updated:
            self._generate_cfg_string()
            self._load_material_data()
            direction = (0,0,1)

        return self._calculate_cross_section(wl, self.matdata, direction=direction)
    
    @staticmethod
    def _get_material_info(material_key: str) -> Dict:
        """Get material information from the materials dictionary."""
        material_info = materials_dict.get(material_key)
        
        if not material_info:
            for info in materials_dict.values():
                if (info.get('formula') == material_key or 
                    info.get('name') == material_key or 
                    info.get('filename') == material_key):
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