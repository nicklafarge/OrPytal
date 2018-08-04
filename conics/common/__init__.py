# Initialize Units
from pint import UnitRegistry, set_application_registry
from conics_utils import orbit_setter
import logging

print('Initializing Unit Registry')
units = UnitRegistry()
Q_ = units.Quantity
set_application_registry(units)

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)