# Initialize Units
from pint import UnitRegistry, set_application_registry

print('Initializing Unit Registry')
ureg = UnitRegistry()
Q_ = ureg.Quantity
set_application_registry(ureg)

