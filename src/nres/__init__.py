"""
nres: Simple yet powerful package for neutron resonance fitting
"""

from __future__ import annotations
from importlib.metadata import version

__all__ = ("__version__",)
__version__ = version(__name__)

from nres.cross_section import CrossSection
from nres.response import Response, Background
from nres.models import TransmissionModel
from nres.data import Data
import nres.utils as utils

materials, elements, isotopes = utils.load_or_create_materials_cache()