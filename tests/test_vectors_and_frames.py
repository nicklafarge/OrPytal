########### Standard ###########
import itertools
import logging
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames, bodies
from orpytal.common import units
from orpytal.planet_constants import CentralBody

########### External ###########
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

earth = bodies.earth

a = 51000 * units.km
i = 1.85 * units.deg
raan = 49.562 * units.deg
argp = 286.537 * units.deg
ta = 23.33 * units.deg

orbit = Orbit(earth, a=a, e=0, raan=raan, arg_periapsis=argp, i=i)
orbit.a = a
orbit.e = 0.2
print(orbit)
print(orbit.angular_momentum)