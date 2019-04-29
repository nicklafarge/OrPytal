########### Standard ###########
import logging

########### Local ###########
from orpytal import Orbit, KeplarianState, get_plot_utils, bodies
from orpytal.common import units

########### External ###########
import numpy as np


logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

# Create an orbit
orbit = Orbit(bodies.earth, e=0.1, p=18400000 * units.m, raan=10*units.deg, inclination=75*units.deg, arg_periapsis=55*units.deg)
state = orbit.get_state(r=18500 * units.km, ascending=True)

# Propagate it
traj = orbit.propagate_orbit()

# Plot it!
oplt = get_plot_utils('plotly')
oplt.init_plot()
oplt.plot_primary(traj)
oplt.plot_traj(traj)
oplt.show()
