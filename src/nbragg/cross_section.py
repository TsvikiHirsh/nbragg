import numpy as np
import pandas as pd
import NCrystal as nc
from typing import Dict, Union, List, Optional
from copy import deepcopy

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

    def __init__(self, materials: Dict[str, Union[Dict, 'CrossSection', float]] = None, 
                 name: str = "", 
                 total_weight: float = 1.,
                 oriented: bool = False):
        """
        Initialize the CrossSection class.

        Args:
            materials: Dictionary of material names and their properties, CrossSection objects, or weights.
            name: Name of the combined cross-section.
            total_weight: Total weight of the cross-section.
            oriented: Whether the materials are oriented or not.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.oriented = oriented

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
                data[column_name] = self.materials[mat](energies) * weight
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
    def from_material(cls, material_file: str, name: str = "", 
                      total_weight: float = 1., temperature: float = 300,
                      oriented: bool = False, mosaicity: float = 0) -> 'CrossSection':
        """
        Create a CrossSection instance from a material file.

        Args:
            material_file: Path to the material file for NCrystal.
            name: Name for the material (optional).
            total_weight: Total weight of the material (default is 1.0).
            temperature: Temperature in Kelvin (default is 300K).
            oriented: Whether the material is oriented or not (default is False).
            mosaicity: Mosaicity in degrees for oriented materials (default is 0).

        Returns:
            CrossSection instance representing the material.
        """
        materials = {
            name: {
                "file": material_file,
                "temperature": temperature,
                "weight": total_weight,
            }
        }
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