name = "orpytal"

from orpytal.common import units, Q_, OrbitType
from orpytal.orbit import Orbit
from orpytal.state import KeplarianState
import orpytal.planet_constants as bodies
import orpytal.ksp.ksp_planet_constants as kerbal_bodies
from orpytal.plotting import get_plot_utils
from orpytal.trajectory import Trajectory

def get_body(body_name):
    return bodies[str(body_name).lower()]