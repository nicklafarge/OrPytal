########### Standard ###########

########### Local ###########
from src.orpytal import Orbit
from common import units

########### External ###########
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit('earth', name='Test Orbit', a=37800*units.km)
print(orbit)