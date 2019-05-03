########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, frames
from orpytal.planet_constants import BODIES_532 as bodies
from orpytal.plotting import plot_orbit
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

earth = bodies['EARTH']
orbit = Orbit(earth, a=10*earth.radius, e=0.4, i=30*units.deg, raan=45*units.deg, arg_periapsis=-90*units.deg)
# st = orbit.get_state(ta=135*units.deg)
st = orbit.get_state()
st.ta=135*units.deg
print(st)