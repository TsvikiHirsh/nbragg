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

import numpy as np
import pandas as pd
import NCrystal as nc

class CrystalCrossSection:
    """
    Represents a combination of cross-sections for crystal materials.
    """

    def __init__(self, materials: Dict[Union[str, 'CrystalCrossSection'], float] = None,
                 name: str = "",
                 total_weight: float = 1.,
                 temp: float = 300.0):
        """
        Initialize the CrystalCrossSection class.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.temp = temp
        self.lambda_grid = np.arange(1.0, 10.0, 0.01)  # Default wavelength grid in Ångstroms
        self.__matdata__ = {}  # Loaded material data
        self.cfg_string = {}  # Configuration strings for each material

        self._generate_cfg_strings()  # Generate configuration strings
        self._load_material_data()  # Load material data
        self._populate_material_data()  # Initialize cross-section data
        self._set_weights()

    def _generate_cfg_strings(self):
        """Generate the configuration strings for all materials."""
        for material in self.materials:
            if isinstance(material, str):
                cfg_string = self._generate_cfg_string(material)
                self.cfg_string[material] = cfg_string

    def _load_material_data(self):
        """Load the material data from .ncmat files using NCrystal."""
        for material, weight in self.materials.items():
            if isinstance(material, str):
                cfg_string = self.cfg_string.get(material)
                if cfg_string:
                    self.__matdata__[material] = nc.createScatter(cfg_string)

    def _generate_cfg_string(self, material, mos=None, temp=None, k=None, l=None):
        """Generate the configuration string for NCrystal based on material type and parameters."""
        mat_info = self._get_material_info(material)
        if not mat_info:
            return ""

        filename = mat_info.get('filename')
        cfg = f"{filename};dcutoff=0.5"

        # Add temperature if provided
        if temp is not None:
            cfg += f";temp={temp}K"

        # Handle single-crystal (oriented) material configuration
        space_group = mat_info.get('space_group', None)
        if space_group is not None and mos is not None:
            cfg += f";mos={mos}deg"
            if k is not None and l is not None:
                cfg += f";dir1=@crys_hkl:{0},{k},{l}@lab:0,0,1;dir2=@crys_hkl:0,-1,1@lab:0,1,0"

        return cfg

    def _get_material_info(self, material_name: str):
        """Retrieve material information from the material dictionary."""
        return materials_dict.get(material_name)

    def _populate_material_data(self):
        """Populate cross-section data for the materials and compute weighted total."""
        xs = {}
        for material, weight in self.materials.items():
            if isinstance(material, str):
                mat_data = self.__matdata__.get(material)
                if mat_data:
                    xs[material] = self._calculate_cross_section(mat_data, self.lambda_grid)
            elif isinstance(material, CrystalCrossSection):
                xs[material.name] = material.table["total"].rename(material.name)

        # Set actual wavelength as index
        self.table = pd.DataFrame(xs, index=self.lambda_grid)
        self.table.index.name = "wavelength"

    def _calculate_cross_section(self, mat_data, lambda_grid):
        """Calculate the cross-section for a material over a wavelength grid."""
        return np.array([mat_data.xsect(wl=λ) for λ in lambda_grid])

    def _set_weights(self):
        """Set and normalize the weights for the materials."""
        self.weights = pd.Series(self.materials)
        self.weights = self.weights[self.weights > 0]
        self.weights /= self.weights.sum()

        self.table["total"] = (self.table * self.weights).sum(axis=1).astype(float)

    def __call__(self, λ: np.ndarray, mos=None, temp=None, k=None, l=None):
        """
        Update material configuration and reload if parameters change.
        Calculate weighted cross-section for a given lambda (λ) array.

        Args:
            λ (np.ndarray): Array of wavelength values (in Ångströms).
        Returns:
            np.ndarray: Weighted cross-section for the given λ values.
        """
        updated = False

        # Check and update material configurations if needed
        for material in self.materials:
            if isinstance(material, str):
                cfg = self.cfg_string.get(material, "")
                new_cfg = self._generate_cfg_string(material, mos, temp, k, l)

                if cfg != new_cfg:
                    self.cfg_string[material] = new_cfg
                    updated = True

        # Reload materials if updated
        if updated:
            self._load_material_data()

        # Calculate and return weighted cross-section for the input λ array
        cross_sections = {}
        for material, weight in self.materials.items():
            mat_data = self.__matdata__.get(material)
            if mat_data:
                cross_sections[material] = np.array([mat_data.xsect(wl=λ_value) for λ_value in λ])

        # Convert to DataFrame for easy manipulation
        cross_section_df = pd.DataFrame(cross_sections, index=λ)
        cross_section_df["total"] = (cross_section_df * self.weights).sum(axis=1)

        # Return weighted total cross-section as np.ndarray
        return cross_section_df["total"].values




    def _lambda_to_energy(self, lambda_value: float) -> float:
        """Convert wavelength in Angstroms to energy in meV."""
        h = 4.135667696e-15  # Planck constant (eV·s)
        m_n = 1.674927471e-27  # Neutron mass (kg)
        energy_meV = (h**2 / (2 * m_n)) * (1 / (lambda_value * 1e-10)**2) * 1e3  # Convert to meV
        return energy_meV

    @classmethod
    def from_material(cls, mat: Union[str, Dict], short_name: str = "", 
                      total_weight: float = 1., temp: float = 300.0) -> 'CrystalCrossSection':
        """
        Create a CrystalCrossSection instance from a material.

        Args:
            mat: Material or dictionary containing material information.
            short_name: Short name for the material (optional).
            total_weight: Total weight of the material (default is 1.0).
            temp: Temperature of the material (in Kelvin).

        Returns:
            CrystalCrossSection instance representing the material.
        """
        if isinstance(mat, str):
            mat_info = cls._get_material_info(mat)
            if not mat_info:
                raise ValueError(f"Material {mat} not found.")

        return cls(materials={mat: total_weight}, name=short_name, temp=temp)

