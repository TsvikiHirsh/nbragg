import numpy as np
import pandas as pd
import NCrystal as nc
from typing import Dict, Union, List, Optional
from copy import deepcopy
import os

from nbragg.utils import materials as materials_dict

class CrossSection:
    """
    Represents a combination of cross-sections for different materials using NCrystal.
    This class loads material cross-section data from NCrystal, allows calculation
    of total cross-section based on weighted sums, and supports both oriented and
    non-oriented materials.

    Attributes:
    materials (Dict[str, Union[Dict, 'CrossSection', float]]): Materials and their properties, CrossSection objects, or weights.
    name (str): Name of the combined cross-section.
    weights (pd.Series): Normalized weights for the materials.
    total_weight (float): Total weight of the cross-section.
    oriented (bool): Whether the materials are oriented or not.
    table (pd.DataFrame): Cross-section data with wavelength as index.
    """

    def __init__(self, materials: Dict[str, Union[Dict, 'CrossSection', float, str]] = None,
                 name: str = "",
                 total_weight: float = 1.,
                 n: float = None,
                 oriented: bool = False):
        """
        Initialize the CrossSection class.

        Args:
        materials: Dictionary of material names and their properties, CrossSection objects, weights, or material names.
        name: Name of the combined cross-section.
        total_weight: Total weight of the cross-section.
        n: atomic density in units of [atoms/barns-cm]
        oriented: Whether the materials are oriented or not.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.oriented = oriented
        self.n = n
        self.nc_materials = {}
        self._load_materials()
        self._set_weights()
        self._generate_table()

    def _load_materials(self):
        """Load the materials using NCrystal or store CrossSection objects."""
        for mat_name, props in self.materials.items():
            if isinstance(props, CrossSection):
                self.nc_materials[mat_name] = props
            elif isinstance(props, (float, int)):
                # If props is a number, assume it's a weight for an existing CrossSection
                self.nc_materials[mat_name] = mat_name  # We'll handle this in __call__
            elif isinstance(props, str):
                # If props is a string, check if it's in materials_dict
                if props in materials_dict:
                    mat_info = materials_dict[props]
                    cfg = f"{mat_info['filename']};temp={300}K"
                    if self.oriented:
                        cfg += f";mos={0}deg"
                    self.nc_materials[mat_name] = nc.load(cfg)
                else:
                    # Assume it's a file path
                    cfg = f"{props};temp={300}K"
                    if self.oriented:
                        cfg += f";mos={0}deg"
                    self.nc_materials[mat_name] = nc.load(cfg)
            elif isinstance(props, dict):
                cfg = f"{props['file']};temp={props.get('temperature', 300)}K"
                if self.oriented:
                    cfg += f";mos={props.get('mosaicity', 0)}deg"
                self.nc_materials[mat_name] = nc.load(cfg)
            else:
                raise ValueError(f"Unsupported material type for {mat_name}")

    def _set_weights(self, weights: Optional[List[float]] = None):
        """
        Set and normalize the weights for the materials.

        Args:
            weights: Optional list of new weights. If provided, must match the number of materials.
        """
        if weights is not None:
            if len(weights) != len(self.materials):
                raise ValueError("Number of weights must match number of materials")
            
            self.weights = pd.Series(weights, index=self.materials.keys())
        else:
            self.weights = pd.Series({mat: props if isinstance(props, (float, int)) else props.get('weight', 1) 
                                      for mat, props in self.materials.items()})

        # Remove materials with zero weight
        self.weights = self.weights[self.weights > 0]

        # Normalize weights
        self.weights /= self.weights.sum()

    def _generate_table(self, wavelengths: np.ndarray = None):
        """
        Generate a pandas DataFrame with cross-section data.
        Args:
        wavelengths: Array of wavelength values. If None, a default range is used.
        """
        if wavelengths is None:
            wavelengths = np.logspace(-2, 1, 1000)  # Default range from 0.01 to 10 Angstroms
        energies = nc.wl2ekin(wavelengths)  # Convert wavelengths to energies
        total_xs = self(energies)
        data = {'total': total_xs}

        for mat_name, weight in self.weights.items():
            mat = self.nc_materials[mat_name]
            column_name = mat_name  # Use the material name as the column name
            if isinstance(mat, CrossSection):
                data[column_name] = mat(energies) * weight
            elif isinstance(mat, str):  # This is the case where we stored the material name
                material = self.materials[mat]
                if callable(material):
                    data[column_name] = material(energies) * weight
                else:
                    # If material is not callable, it's likely a weight or other scalar
                    data[column_name] = np.full_like(energies, material * weight)
            else:
                if self.oriented:
                    absorption = mat.absorption.crossSection(energies)
                    scatter = mat.scatter.crossSection(energies)
                else:
                    absorption = mat.absorption.crossSectionNonOriented(energies)
                    scatter = mat.scatter.crossSectionNonOriented(energies)
                data[column_name] = (absorption + scatter) * weight

        self.table = pd.DataFrame(data, index=wavelengths)
        self.table.index.name = 'wavelength'



    def __call__(self, E: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Calculate the weighted cross-section for a given set of energies.

        Args:
            E: Array of energy values.
            weights: Optional array of new weights.

        Returns:
            Array of weighted cross-section values.
        """
        if weights is not None:
            self._set_weights(weights=weights)

        total_xs = np.zeros_like(E)
        for mat_name, weight in self.weights.items():
            mat = self.nc_materials[mat_name]
            if isinstance(mat, CrossSection):
                total_xs += weight * mat(E)
            elif isinstance(mat, str):  # This is the case where we stored the material name
                total_xs += weight * self.materials[mat](E)
            else:
                if self.oriented:
                    absorption = mat.absorption.crossSection(E)
                    scatter = mat.scatter.crossSection(E)
                else:
                    absorption = mat.absorption.crossSectionNonOriented(E)
                    scatter = mat.scatter.crossSectionNonOriented(E)
                total_xs += weight * (absorption + scatter)

        return total_xs
    
    def __add__(self, other: 'CrossSection') -> 'CrossSection':
        """
        Add two CrossSection objects.
        Args:
        other: Another CrossSection object to add to the current one.
        Returns:
        A new CrossSection object representing the sum of the two.
        """
        if self.oriented != other.oriented:
            raise ValueError("Cannot add oriented and non-oriented CrossSections")

        new_self = CrossSection()  # Create a new CrossSection object
        
        # Combine materials
        new_self.materials = {**self.materials, **other.materials}
        
        # Combine weights
        self_weights = pd.Series(self.weights, name='weight')
        other_weights = pd.Series(other.weights, name='weight')
        new_weights = (self_weights * self.total_weight).add(
            other_weights * other.total_weight, fill_value=0
        )
        new_self.weights = new_weights / new_weights.sum()
        
        new_self.total_weight = 1.
        new_self.oriented = self.oriented
        
        # Copy other necessary attributes (add any that are missing)
        new_self.nc_materials = {**self.nc_materials, **other.nc_materials}
        
        new_self._load_materials()
        new_self._generate_table()
        return new_self

    def plot(self, wavelengths: np.ndarray = None, **kwargs):
        """
        Plot the cross-section data.

        Args:
            wavelengths: Array of wavelength values. If None, uses the table's index.
            **kwargs: Optional plotting parameters.
        """
        import matplotlib.pyplot as plt

        if wavelengths is not None:
            self._generate_table(wavelengths)

        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "$\sigma$ [barn]")
        xlabel = kwargs.pop("xlabel", "Wavelength [Å]")
        lw = kwargs.pop("lw", 1.)

        fig, ax = plt.subplots()
        self.table.plot(ax=ax, logy=True, logx=True)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        return ax

    @classmethod
    def from_material(cls, material: Union[str, Dict], name: str = "",
                      total_weight: float = 1., temperature: float = 300,
                      oriented: bool = False, mosaicity: float = 0,
                      atomic_density: float = None) -> 'CrossSection':
        """
        Create a CrossSection instance from a material file, material dictionary, or material name.

        Args:
        material: Path to the material file for NCrystal, a dictionary containing material information, or a material name.
        name: Name for the material (optional).
        total_weight: Total weight of the material (default is 1.0).
        temperature: Temperature in Kelvin (default is 300K).
        oriented: Whether the material is oriented or not (default is False).
        mosaicity: Mosaicity in degrees for oriented materials (default is 0).
        atomic_density: Atomic density of the material in units of [atoms/barn-cm]

        Returns:
        CrossSection instance representing the material.
        """
        if isinstance(material, str):
            # If material is a string, check if it's in materials_dict
            if material in materials_dict:
                material = materials_dict[material]
                name = name or material['name']
                materials = {
                    name: {
                        "file": material['filename'],
                        "temperature": temperature,
                        "weight": total_weight,
                        "n": material.get("atomic_density", atomic_density)
                    }
                }
            else:
                # Assume it's a file path
                materials = {
                    name: {
                        "file": material,
                        "temperature": temperature,
                        "weight": total_weight,
                        "n": atomic_density
                    }
                }
        elif isinstance(material, dict):
            # If material is a dictionary, use it directly
            name = material.get('name', name)
            materials = {
                name: {
                    "file": material['filename'],
                    "temperature": temperature,
                    "weight": total_weight,
                    "n": material.get("atomic_density", atomic_density)
                }
            }
        else:
            raise ValueError("material must be either a file path string, a material name, or a material dictionary")

        if oriented:
            materials[name]["mosaicity"] = mosaicity

        return cls(materials, name=name, total_weight=total_weight, oriented=oriented)
    
    def update_weights(self, new_weights: Dict[str, float]):
        """
        Update the weights of the materials.

        Args:
            new_weights: Dictionary of material names and their new weights.
        """
        for mat, weight in new_weights.items():
            if mat not in self.weights.index:
                raise ValueError(f"Material {mat} not found in the CrossSection")
            self.weights[mat] = weight
        
        self._set_weights()
        self._generate_table()

    def add_material(self, material: Union[Dict, 'CrossSection'], name: str, weight: float):
        """
        Add a new material to the CrossSection.

        Args:
            material: Material properties dictionary or CrossSection object.
            name: Name of the new material.
            weight: Weight of the new material.
        """
        self.materials[name] = material
        self._load_materials()
        self.weights[name] = weight
        self._set_weights()
        self._generate_table()

    def remove_material(self, name: str):
        """
        Remove a material from the CrossSection.

        Args:
            name: Name of the material to remove.
        """
        if name not in self.materials:
            raise ValueError(f"Material {name} not found in the CrossSection")
        
        del self.materials[name]
        del self.nc_materials[name]
        self.weights = self.weights.drop(name)
        self._set_weights()
        self._generate_table()

    def get_material_xs(self, material_name: str, wavelengths: np.ndarray = None) -> np.ndarray:
        """
        Get the cross-section for a specific material.

        Args:
            material_name: Name of the material.
            wavelengths: Array of wavelength values. If None, uses the table's index.

        Returns:
            Array of cross-section values for the specified material.
        """
        if material_name not in self.materials:
            raise ValueError(f"Material {material_name} not found in the CrossSection")

        if wavelengths is not None:
            self._generate_table(wavelengths)

        return self.table[material_name].values

    def get_total_xs(self, wavelengths: np.ndarray = None) -> np.ndarray:
        """
        Get the total cross-section.

        Args:
            wavelengths: Array of wavelength values. If None, uses the table's index.

        Returns:
            Array of total cross-section values.
        """
        if wavelengths is not None:
            self._generate_table(wavelengths)

        return self.table['total'].values

import os
import pandas as pd
import numpy as np
import NCrystal as nc
from typing import Dict, Union, Optional
from copy import deepcopy

class CrystalCrossSection:
    """
    Represents a combination of cross-sections for crystal materials.
    """

    def __init__(self, materials: Dict[Union[str, 'CrystalCrossSection'], float] = None,
                name: str = "",
                total_weight: float = 1.0,
                temp: float = 300.0,
                mos: Union[float, Dict[str, float], None] = None,
                k: Union[float, Dict[str, float], None] = None,
                l: Union[float, Dict[str, float], None] = None):
        """
        Initialize the CrystalCrossSection class.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.temp = temp
        self.mos = self._process_parameter(mos)
        self.k = self._process_parameter(k)
        self.l = self._process_parameter(l)
        
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

            cfg += f";mos={mos}deg"
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
            elif isinstance(material, CrystalCrossSection):
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
                      mos=None, k=None, l=None) -> 'CrystalCrossSection':
        """
        Create a CrystalCrossSection instance from a single material.

        Args:
            mat: Material name or dictionary with material info.
            short_name: Optional short name for the material.
            total_weight: Total weight for the material.
            temp: Temperature (in K).
            mos: Mosaicity of the material.
            k, l: Orientation indices.
        
        Returns:
            A CrystalCrossSection instance for the material.
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

        materials = {filename: total_weight}
        return cls(materials=materials, name=short_name, total_weight=total_weight, 
                   temp=temp, mos=mos, k=k, l=l)

    def __add__(self, other: 'CrystalCrossSection') -> 'CrystalCrossSection':
        """
        Add two CrystalCrossSection objects, combining their materials and weights.

        Args:
            other: Another CrystalCrossSection object.
        
        Returns:
            A new CrystalCrossSection object representing the sum.
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
        return CrystalCrossSection(materials=new_weights_normalized.to_dict(), 
                                total_weight=new_weights.sum(), 
                                mos=new_mos, k=new_k, l=new_l, 
                                temp=self.temp)

    def __mul__(self, scalar: float) -> 'CrystalCrossSection':
        """
        Multiply the total weight of the CrystalCrossSection by a scalar.

        Args:
            scalar: The value to multiply the total weight by.

        Returns:
            A new CrystalCrossSection object with the updated material weights and total weight.
        """
        # Scale the material weights by the scalar
        scaled_weights = pd.Series(self.materials) * scalar
        
        # No need to normalize here; the user may want to combine with other objects later.
        return CrystalCrossSection(materials=scaled_weights.to_dict(), 
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