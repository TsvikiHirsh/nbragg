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
                dirtol:float = 1.):
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
        
        self.lambda_grid = np.arange(1.0, 10.0, 0.01)  # Default wavelength grid in Ångstroms
        self.__matdata__ = {}  # Loaded material data
        self.cfg_string = {}  # Configuration strings for each material

        self._set_weights()
        self._generate_cfg_strings()  # Generate configuration strings
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

    def _generate_cfg_strings(self):
        for material in self.materials:
            if isinstance(material, str):
                cfg_string = self._generate_cfg_string(material)
                self.cfg_string[material] = cfg_string

    def _generate_cfg_string(self, material):
        mat_info = self._get_material_info(material)
        if not mat_info:
            return ""

        filename = mat_info.get('filename')
        cfg = f"{filename};dcutoff=0.5"

        if self.temp is not None:
            cfg += f";temp={self.temp}K"

        mos = self.mos.get(material)
        k = self.k.get(material)
        l = self.l.get(material)

        if any(param is not None for param in (mos, k, l)):
            mos = mos if mos is not None else 0.001
            k = k if k is not None else 0
            l = l if l is not None else 1

            cfg += f";mos={mos}deg;dirtol={self.dirtol}deg"
            cfg += f";dir1=@crys_hkl:0,{k},{l}@lab:0,0,1"
            cfg += f";dir2=@crys_hkl:0,-1,1@lab:0,1,0"

        return cfg

    def _load_material_data(self):
        """Load the material data from .ncmat files using NCrystal."""
        for material, weight in self.materials.items():
            if isinstance(material, str):
                cfg_string = self.cfg_string.get(material)
                if cfg_string:
                    self.__matdata__[material] = nc.createScatter(cfg_string)

    def _get_material_info(self, material_name: str):
        """Retrieve material information from the material dictionary."""
        return materials_dict.get(material_name)

    def _set_weights(self):
        self.weights = pd.Series(self.materials)
        self.weights = self.weights[self.weights > 0]
        self.weights /= self.weights.sum()

    def _populate_material_data(self):
        xs = {}
        for material, weight in self.materials.items():
            if isinstance(material, str):
                mat_data = self.__matdata__.get(material)
                if mat_data:
                    xs[material] = self._calculate_cross_section(self.lambda_grid, mat_data, material)
            elif isinstance(material, CrossSection):
                xs[material.name] = material.table["total"].rename(material.name)

        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"
        self.table["total"] = (self.table * self.weights).sum(axis=1).astype(float)


    def _calculate_cross_section(self, λ, mat_data, material):
        ekin = nc.wl2ekin(λ)
        mos = self.mos.get(material)
        k = self.k.get(material)
        l = self.l.get(material)

        if all(param is not None for param in (mos, k, l)):
            # Use the oriented cross-section
            return np.array([
                mat_data.crossSection(ekin_val, direction=(0, 0, 1))
                for ekin_val in ekin
            ])
        else:
            # Use the non-oriented cross-section
            return mat_data.crossSectionNonOriented(ekin)
        

    def __call__(self, λ: np.ndarray, mos=None, temp=None, k=None, l=None):
        """
        Update material configuration if parameters change and return fast calculation of
        the weighted cross-section for a given lambda (λ) array.
        
        Args:
            λ (np.ndarray): Array of wavelength values (in Ångströms).
            
        Returns:
            np.ndarray: Weighted cross-section for the given λ values.
        """
        updated = False

        # Update parameters if provided
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

        # Regenerate cfg strings and reload material data if updated
        if updated:
            self._generate_cfg_strings()
            self._load_material_data()

        # Calculate and return weighted cross-section for the input λ array
        cross_sections = {}
        for material, weight in self.materials.items():
            mat_data = self.__matdata__.get(material)
            if mat_data:
                cross_sections[material] = self._calculate_cross_section(λ, mat_data, material)

        # Convert to DataFrame for easy manipulation
        cross_section_df = pd.DataFrame(cross_sections, index=λ)
        
        # Multiply by weights and sum the cross-sections for all materials
        total_cross_section = (cross_section_df * self.weights).sum(axis=1)

        # Return the weighted total cross-section as np.ndarray
        return total_cross_section.values



    def _lambda_to_energy(self, lambda_value: float) -> float:
        """Convert wavelength in Angstroms to energy in meV."""
        h = 4.135667696e-15  # Planck constant (eV·s)
        m_n = 1.674927471e-27  # Neutron mass (kg)
        energy_meV = (h**2 / (2 * m_n)) * (1 / (lambda_value * 1e-10)**2) * 1e3  # Convert to meV
        return energy_meV

    @classmethod
    def from_material(cls, mat: Union[str, Dict], short_name: str = "", 
                      total_weight: float = 1.0, temp: float = 300.0, 
                      mos=None, k=None, l=None,dirtol=None) -> 'CrossSection':
        """
        Create a CrossSection instance from a single material.

        Args:
            mat: Material name or dictionary with material info.
            short_name: Optional short name for the material.
            total_weight: Total weight for the material.
            temp: Temperature (in K).
            mos: Mosaicity of the material.
            k, l: Orientation indices.
        
        Returns:
            A CrossSection instance for the material.
        """
        if isinstance(mat, str):
            mat_info = cls._get_material_info(mat)
            if not mat_info:
                raise ValueError(f"Material '{mat}' not found.")
        elif isinstance(mat, dict):
            mat_info = mat
        else:
            raise TypeError("Argument 'mat' must be a string or dictionary.")

        short_name = mat_info.get('short_name', short_name)
        filename = mat_info.get('filename')
        if not filename:
            raise ValueError(f"Material '{mat}' does not contain a valid filename.")
        
        dirtol = dirtol if dirtol else 1.

        materials = {filename: total_weight}
        return cls(materials=materials, name=short_name, total_weight=total_weight, 
                   temp=temp, mos=mos, k=k, l=l, dirtol=dirtol)

    def __add__(self, other: 'CrossSection') -> 'CrossSection':
        """
        Add two CrossSection objects, combining their materials and weights.

        Args:
            other: Another CrossSection object.
        
        Returns:
            A new CrossSection object representing the sum.
        """
        # Combine materials
        new_materials = {**self.materials, **other.materials}
        
        # Get the weighted materials for self and other
        self_weights = pd.Series(self.materials) * self.total_weight
        other_weights = pd.Series(other.materials) * other.total_weight
        
        # Add the weights and normalize them to sum to 1
        new_weights = self_weights.add(other_weights, fill_value=0)
        new_weights_normalized = new_weights / new_weights.sum()
        
        # Combine mos, k, l values
        new_mos = self._combine_dict_or_list(self.mos, other.mos)
        new_k = self._combine_dict_or_list(self.k, other.k)
        new_l = self._combine_dict_or_list(self.l, other.l)

        # Create a new object with the combined attributes
        return CrossSection(materials=new_weights_normalized.to_dict(), 
                                total_weight=new_weights.sum(), 
                                mos=new_mos, k=new_k, l=new_l, 
                                temp=self.temp)

    def __mul__(self, scalar: float) -> 'CrossSection':
        """
        Multiply the total weight of the CrossSection by a scalar.

        Args:
            scalar: The value to multiply the total weight by.

        Returns:
            A new CrossSection object with the updated material weights and total weight.
        """
        # Scale the material weights by the scalar
        scaled_weights = pd.Series(self.materials) * scalar
        
        # No need to normalize here; the user may want to combine with other objects later.
        return CrossSection(materials=scaled_weights.to_dict(), 
                                total_weight=self.total_weight * scalar, 
                                mos=self.mos, k=self.k, l=self.l, 
                                temp=self.temp)


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
                if info['formula'] == material_key or info['short_name'] == material_key or info['filename'] == material_key:
                    return info

        return material_info
    
    def plot(self, **kwargs):
        """
        Plot the cross-section data.

        Args:
            wavelengths: Array of wavelength values. If None, uses the table's index.
            **kwargs: Optional plotting parameters.
        """
        import matplotlib.pyplot as plt

        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "$\sigma$ [barn]")
        xlabel = kwargs.pop("xlabel", "Wavelength [Å]")
        lw = kwargs.pop("lw", 1.)

        fig, ax = plt.subplots()
        self.table.mul(np.r_[self.weights,1]).plot(ax=ax, logy=True, logx=True)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        return ax