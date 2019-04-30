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

orbit = Orbit(bodies['EARTH'], se=0, rp=6603)
st = orbit.get_state(r=2*orbit.central_body.radius)
print(st)