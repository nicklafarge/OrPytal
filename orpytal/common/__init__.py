
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
    # Circular = 0
    Elliptic = 0
    Parabolic = 1
    Hyperbolic = 2