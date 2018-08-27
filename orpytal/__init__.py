name = "orpytal"

from orpytal.orbit import Orbit
from orpytal.state import KeplarianState
from orpytal.common import units, Q_
import orpytal.planet_constants as bodies

def get_body(body_name):
    return bodies[str(body_name).lower()]