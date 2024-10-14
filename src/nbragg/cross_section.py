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
        materials (Dict[str, Dict]): Materials and their properties.
        name (str): Name of the combined cross-section.
        weights (pd.Series): Normalized weights for the materials.
        total_weight (float): Total weight of the cross-section.
        oriented (bool): Whether the materials are oriented or not.
    """

    def __init__(self, materials: Dict[str, Dict] = None, 
                 name: str = "", 
                 total_weight: float = 1.,
                 oriented: bool = False):
        """
        Initialize the CrossSection class.

        Args:
            materials: Dictionary of material names and their properties.
            name: Name of the combined cross-section.
            total_weight: Total weight of the cross-section.
            oriented: Whether the materials are oriented or not.
        """
        self.materials = materials or {}
        self.name = name
        self.total_weight = total_weight if self.materials else 0.
        self.oriented = oriented

        self._load_materials()
        self._set_weights()

    def _load_materials(self):
        """Load the materials using NCrystal."""
        self.nc_materials = {}
        for mat_name, props in self.materials.items():
            cfg = f"{props['file']};temp={props.get('temperature', 300)}K"
            if self.oriented:
                cfg += f";mos={props.get('mosaicity', 0)}deg"
            self.nc_materials[mat_name] = nc.load(cfg)

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
            self.weights = pd.Series({mat: props.get('weight', 1) for mat, props in self.materials.items()})

        # Remove materials with zero weight
        self.weights = self.weights[self.weights > 0]

        # Normalize weights
        self.weights /= self.weights.sum()

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

        new_materials = {**self.materials, **other.materials}
        new_weights = (self.weights * self.total_weight).add(
            other.weights * other.total_weight, fill_value=0
        )
        new_weights /= new_weights.sum()

        new_self = deepcopy(self)
        new_self.materials = new_materials
        new_self.weights = new_weights
        new_self.total_weight = 1.
        new_self._load_materials()

        return new_self

    def __mul__(self, total_weight: float = 1.) -> 'CrossSection':
        """
        Multiply the CrossSection by a total weight.

        Args:
            total_weight: The weight to multiply by.

        Returns:
            A new CrossSection object with updated total_weight.
        """
        new_self = deepcopy(self)
        new_self.total_weight = total_weight
        return new_self

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
            if self.oriented:
                absorption = mat.absorption.crossSection(E)
                scatter = mat.scatter.crossSection(E)
            else:
                absorption = mat.absorption.crossSectionNonOriented(E)
                scatter = mat.scatter.crossSectionNonOriented(E)
            total_xs += weight * (absorption + scatter)

        return total_xs

    def plot(self, E: np.ndarray, **kwargs):
        """
        Plot the cross-section data.

        Args:
            E: Array of energy values.
            **kwargs: Optional plotting parameters.
        """
        import matplotlib.pyplot as plt

        title = kwargs.pop("title", self.name)
        ylabel = kwargs.pop("ylabel", "$\sigma$ [barn]")
        xlabel = kwargs.pop("xlabel", "Energy [eV]")
        lw = kwargs.pop("lw", 1.)

        fig, ax = plt.subplots()
        total_xs = self(E)
        ax.plot(E, total_xs, label="Total", linewidth=1.5, color="0.2", zorder=100)

        for mat_name, weight in self.weights.items():
            mat = self.nc_materials[mat_name]
            if self.oriented:
                absorption = mat.absorption.crossSection(E)
                scatter = mat.scatter.crossSection(E)
            else:
                absorption = mat.absorption.crossSectionNonOriented(E)
                scatter = mat.scatter.crossSectionNonOriented(E)
            xs = weight * (absorption + scatter)
            ax.plot(E, xs, label=f"{mat_name}: {weight*100:>6.2f}%", linewidth=lw)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.set_xscale('log')
        ax.set_yscale('log')
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