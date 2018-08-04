########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState, plotting
from common import units, Q_

########### External ###########
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit('earth', name='Test Orbit', a=37800*units.km)
print(orbit)