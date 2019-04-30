########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils
from orpytal.planet_constants import BODIES_532 as bodies
from orpytal.plotting import plot_orbit
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

# earth = bodies['EARTH']
orbit = Orbit(bodies['MOON'], rp=1937.4, a=-6460)
st = orbit.get_state(ta=-90*units.deg)
print(st)

print(orbit.get_state(ta=0).v)