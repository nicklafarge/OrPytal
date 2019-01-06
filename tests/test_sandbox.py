########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils

import matplotlib.pyplot as plt
import numpy as np

# logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)


orbit = Orbit(bodies.earth, a=51000 * units.km, e=0.7, raan=10*units.deg, arg_periapsis=10*units.deg, inclination=45*units.deg)
st = KeplarianState(orbit)
# st.r = 10 * units.deg
# st.v = 1 * units("km/s")
st.t_since_rp = 2 * units.day
print(st)

# traj = orbit.propagate_full_orbit()
# pu = get_plot_utils("plotly", planar=False)
# pu.init_plot()
# pu.plot_traj(traj)
# pu.plot_primary(bodies.earth)
# pu.show()