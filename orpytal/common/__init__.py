# Initialize Units
from pint import UnitRegistry, set_application_registry

units = UnitRegistry()
Q_ = units.Quantity
set_application_registry(units)

# Initialize Logging
import logging

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

from enum import Enum


class OrbitType(Enum):
    Unknown = 0
    Elliptic = 1
    Parabolic = 2
    Hyperbolic = 3
