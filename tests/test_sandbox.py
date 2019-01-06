########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils

import matplotlib.pyplot as plt
import numpy as np

# logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)

orbit = Orbit(bodies.earth, a=51000, e=0.0)
pair = ('a', 'e')

test_orbit = Orbit(orbit.central_body)
setattr(test_orbit, pair[0], getattr(orbit, pair[0]))
setattr(test_orbit, pair[1], getattr(orbit, pair[1]))
same = orbit.compare(test_orbit)