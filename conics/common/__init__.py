# Initialize Units
from pint import UnitRegistry, set_application_registry
from conics_utils import orbit_setter

print('Initializing Unit Registry')
ureg = UnitRegistry()
Q_ = ureg.Quantity
set_application_registry(ureg)

