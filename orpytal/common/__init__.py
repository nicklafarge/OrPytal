# Initialize Units
from pint import UnitRegistry, set_application_registry
from orpytal.common.conics_utils import orbit_setter, attribute_setter

import logging

units = UnitRegistry()
Q_ = units.Quantity
set_application_registry(units)

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)