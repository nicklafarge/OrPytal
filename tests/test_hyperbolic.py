########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils

import matplotlib.pyplot as plt
import numpy as np

# logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)

orbit = Orbit(bodies.earth, a=-51000, e=1.1)
state = KeplarianState(orbit, ta=30 * units.deg)
print(state)