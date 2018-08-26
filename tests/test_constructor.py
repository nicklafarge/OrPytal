########### Standard ###########

########### Local ###########
from orpytal import Orbit, bodies, plotting
from orpytal.common import units

########### External ###########
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit(bodies.earth, name='Test Orbit', a=37800*units.km, p=16000*units.km)
plotting.plot_orbit(orbit)

print(orbit)