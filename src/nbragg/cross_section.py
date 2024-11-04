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
    def __init__(self, materials: Dict[Union[str, 'CrossSection'], float] = None,
                name: str = "",
                total_weight: float = 1.0,
                temp: float = 300.0,
                mos: Union[float, Dict[str, float], None] = None,
                k: Union[float, Dict[str, float], None] = None,
                l: Union[float, Dict[str, float], None] = None,
                dirtol: float = 1.):
        """
        Initialize the CrossSection class.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.temp = temp
        self.mos = self._process_parameter(mos)
        self.k = self._process_parameter(k)
        self.l = self._process_parameter(l)
        self.dirtol = dirtol
        self.L = 9.  # TODO: replace this hack
        
        self.lambda_grid = np.arange(1.0, 10.0, 0.01)  # Default wavelength grid in Ångstroms
        self.matdata = None  # Single NCrystal scatter object
        
        # Initialize weights as empty Series first
        self.weights = pd.Series(dtype=float)
        self._set_weights()
        self._generate_cfg_string()  # Generate single configuration string
        self._load_material_data()  # Load material data
        self._populate_material_data()  # Initialize cross-section data
        

    def _process_parameter(self, param):
        if isinstance(param, (int, float)):
            return {m: param for m in self.materials}
        elif isinstance(param, (list, np.ndarray)):
            return {m: param[i] if i < len(param) else None for i, m in enumerate(self.materials)}
        elif isinstance(param, dict):
            return param
        else:
            return {m: None for m in self.materials}

    def _generate_cfg_string(self):
        """Generate a single configuration string using NCrystal phase notation."""
        if self.weights is None or self.weights.empty:
            self.cfg_string = ""
            return
            
        # Build the phase string
        phase_parts = []
        
        for material, weight in self.weights.items():
            if isinstance(material, str):
                mat_info = self._get_material_info(material)
                if not mat_info:
                    continue
                    
                filename = mat_info.get('filename')
                phase = f"{weight}*{filename}"
                
                # Add material-specific parameters
                params = []
                
                if self.temp is not None:
                    params.append(f"temp={self.temp}K")
                
                mos = self.mos.get(material)
                k = self.k.get(material)
                l = self.l.get(material)
                
                if all(param is not None for param in (mos, k, l)):
                    params.append(f"mos={mos}deg")
                    params.append(f"dirtol={self.dirtol}deg")
                    params.append(f"dir1=@crys_hkl:0,{k},{l}@lab:0,0,1")
                    params.append(f"dir2=@crys_hkl:0,-1,1@lab:0,1,0")
                
                if params:
                    phase += f";{';'.join(params)}"
                
                phase_parts.append(phase)
            
            elif isinstance(material, CrossSection):
                # For nested CrossSection objects, include their cfg_string with proper weight
                if hasattr(material, 'cfg_string') and material.cfg_string:
                    phase_parts.append(f"{weight}*({material.cfg_string})")
        
        # Combine all phases
        self.cfg_string = f"phases<{'&'.join(phase_parts)}>" if phase_parts else ""
        self.phases = phase_parts



    def _load_material_data(self):
        """Load the material data using NCrystal with the phase configuration."""
        if self.cfg_string:
            self.matdata = nc.load(self.cfg_string)
            self.phases = {name:{"mat":nc.NCMATComposer.from_info(phase[1]).load(),"weight":phase[0]} for name,phase in zip(self.phases,self.matdata.info.phases)}


    def _get_material_info(self, material_name: str):
        """Retrieve material information from the material dictionary."""
        return materials_dict.get(material_name)
    
    def _set_weights(self, name: str = ""):
        """Set weights for materials and ensure they're properly normalized."""
        name = name if name else self.name
        materials_with_names = {}
        
        for key, value in self.materials.items():
            if isinstance(key, str):
                # For string keys, use the material name or the provided name
                materials_with_names[name if key == self.materials.get(key) else key] = value
            elif isinstance(key, CrossSection):
                # For CrossSection objects, use their name
                materials_with_names[key.name] = value
            else:
                # For other cases, use the key directly
                materials_with_names[key] = value
        
        # Create a pandas Series with the weights
        if materials_with_names:
            weights = pd.Series(materials_with_names)
            # Filter out zero or negative weights
            weights = weights[weights > 0]
            # Normalize weights to sum to 1
            if not weights.empty:
                weights = weights / weights.sum()# * self.total_weight
            self.weights = weights
        else:
            # Initialize empty weights if no materials
            self.weights = pd.Series(dtype=float)



    def _populate_material_data(self):
        """Populate cross section data using NCrystal phases."""
        if not self.cfg_string:
            self.table = pd.DataFrame(index=self.lambda_grid)
            self.table.index.name = "wavelength"
            return

        # Load the material with all phases
        mat = nc.load(self.cfg_string)
        xs = {}

        # Process each phase separately
        for phase, phase_info in self.phases.items():
            # Calculate cross-section for this phase
            xs[phase] = self._calculate_cross_section(self.lambda_grid, phase_info["mat"])

        # calculate total
        xs["total"] = self._calculate_cross_section(self.lambda_grid, self.matdata) 

        # Create DataFrame with all phases
        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"
        if len(self.table.columns)>1:
            self.table.columns = self.weights.index.to_list() + ["total"]
        else:
            self.table.columns = ["total"]

        self.atomic_density = mat.info.factor_macroscopic_xs


    def _calculate_cross_section(self, wl, mat, direction = None):
        """
        Calculate cross-section using NCrystal's xsect method.
        
        Args:
            wl (array-like): Wavelength values in Angstroms
            material (str): Material name for parameter lookup
            
        Returns:
            numpy.ndarray: Array of cross-section values
        """
        return mat.scatter.xsect(wl=wl, direction=direction) + mat.absorption.xsect(wl=wl, direction=direction)



        


    def __call__(self, wl: np.ndarray, mos=None, temp=None, k=None, l=None):
        """Update configuration if parameters change and return cross-section."""
        updated = False
        direction = None
        if mos is not None:
            self.mos = self._process_parameter(mos)
            updated = True
        if temp is not None:
            self.temp = temp
            updated = True
        if k is not None:
            self.k = self._process_parameter(k)
            updated = True
        if l is not None:
            self.l = self._process_parameter(l)
            updated = True

        if updated:
            self._generate_cfg_string()
            self._load_material_data()
            direction = (0,0,1)

        return self._calculate_cross_section(wl,self.matdata,direction=direction)



    def _wavelength_to_energy(self, wavelength: float) -> float:
        """Convert wavelength in Angstroms to energy in meV."""
        h = 4.135667696e-15  # Planck constant (eV·s)
        m_n = 1.674927471e-27  # Neutron mass (kg)
        energy_meV = (h**2 / (2 * m_n)) * (1 / (wavelength * 1e-10)**2) * 1e3  # Convert to meV
        return energy_meV

    @classmethod
    def from_material(cls, mat: Union[str, Dict], short_name: str = "", 
                      total_weight: float = 1.0, temp: float = 300.0, 
                      mos=None, k=None, l=None, dirtol=None) -> 'CrossSection':
        """
        Create a CrossSection instance from a single material.
        """
        if isinstance(mat, str):
            mat_info = materials_dict.get(mat)
            if not mat_info:
                raise ValueError(f"Material '{mat}' not found.")
        elif isinstance(mat, dict):
            mat_info = mat
        else:
            raise TypeError("Argument 'mat' must be a string or dictionary.")

        # Use the material's short name if none provided
        short_name = short_name if short_name != "" else mat_info.get('short_name', '')
        
        # Get the filename from material info
        filename = mat_info.get('filename')
        if not filename:
            raise ValueError(f"Material '{mat}' does not contain a valid filename.")
        
        dirtol = dirtol if dirtol is not None else 1.

        # Create materials dictionary with the filename as key
        materials = {filename: total_weight}
        
        # Create instance
        instance = cls(materials=materials, 
                      name=short_name, 
                      total_weight=total_weight,
                      temp=temp, 
                      mos=mos, 
                      k=k, 
                      l=l, 
                      dirtol=dirtol)
        
        # Update weights with proper index
        if short_name:
            instance.weights.index = [short_name]
        
        return instance
    
    def __add__(self, other: 'CrossSection') -> 'CrossSection':
        """Add two CrossSection objects."""
        new_materials = {**self.materials, **other.materials}
        new_weights = {}
        
        # Combine weights
        for material, weight in self.materials.items():
            new_weights[material] = weight * self.total_weight
        for material, weight in other.materials.items():
            if material in new_weights:
                new_weights[material] += weight * other.total_weight
            else:
                new_weights[material] = weight * other.total_weight
        
        # Combine parameters
        new_mos = self._combine_dict_or_list(self.mos, other.mos)
        new_k = self._combine_dict_or_list(self.k, other.k)
        new_l = self._combine_dict_or_list(self.l, other.l)
        
        instance = CrossSection(
            materials=new_weights,
            total_weight=sum(new_weights.values()),
            mos=new_mos,
            k=new_k,
            l=new_l,
            temp=self.temp
        )
        
        return instance

    def __mul__(self, scalar: float) -> 'CrossSection':
        """Multiply CrossSection by a scalar."""
        new_materials = {m: w * scalar for m, w in self.materials.items()}
        new_total_weight = self.total_weight * scalar
        
        # Normalize the new weights
        new_weights = pd.Series(new_materials)
        new_weights = new_weights[new_weights > 0]
        if not new_weights.empty:
            new_weights = new_weights / new_weights.sum()
        
        return CrossSection(
            materials=new_weights.to_dict(),
            total_weight=new_total_weight,
            mos=self.mos,
            k=self.k,
            l=self.l,
            temp=self.temp
        )


    def _combine_dict_or_list(self, a, b):
        """
        Helper function to combine dictionaries or lists.
        If both are dictionaries, keys will be merged.
        If both are lists, values will be concatenated.
        """
        if isinstance(a, dict) and isinstance(b, dict):
            return {**a, **b}
        if isinstance(a, list) and isinstance(b, list):
            return a + b
        return a  # Default behavior if they are not both dicts or lists

    @classmethod
    def _get_material_info(cls, material_key: str) -> Dict:
        """
        Get material information by looking up the material by formula, short_name, or filename.

        Args:
            material_key: Either the formula, short name, or filename of the material.

        Returns:
            A dictionary containing the material information if found, otherwise None.
        """
        # Search by formula, short name, or filename in materials_dict
        material_info = materials_dict.get(material_key)

        if not material_info:
            # If not found directly, search by matching values in the materials dictionary
            for key, info in materials_dict.items():
                if info['formula'] == material_key or info['name'] == material_key or info['filename'] == material_key:
                    return info

        return material_info
    
    def plot(self, **kwargs):
        """
        Plot the cross-section data.

        Args:
            **kwargs: Optional plotting parameters.
        """
        import matplotlib.pyplot as plt
        
        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "$\sigma$ [barn]")
        xlabel = kwargs.pop("xlabel", "Wavelength [Å]")
        lw = kwargs.pop("lw", 1.)

        fig, ax = plt.subplots()

        # Plot each material component with reduced line width
        if len(self.table.columns)>1:
            self.table.iloc[:, :-1].mul(self.weights).plot(ax=ax, lw=lw, **kwargs)

        # Plot the total curve with a thicker line width and distinct color
        self.table["total"].plot(ax=ax, color="0.2", lw=lw*1.2, label="Total")

        # Set labels, title, and legend
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        # Legend with formatted names and weights for components plus 'total'
        legend_labels = [f"{material}: {weight*100:.1f}%" for material, weight in self.weights.items()] + ["Total"]
        ax.legend(legend_labels)

        return ax