"""
nbragg: Simple yet powerful package for neutron resonance fitting
"""

from __future__ import annotations
from importlib.metadata import version

__all__ = ("__version__",)
__version__ = version(__name__)

from nbragg.cross_section import CrossSection
from nbragg.response import Response, Background
from nbragg.models import TransmissionModel
from nbragg.data import Data
import nbragg.utils as utils
from nbragg.utils import materials, register_material


