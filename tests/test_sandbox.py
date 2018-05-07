########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState, plotting
from common import ureg, Q_

########### External ###########
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit('earth', name='Test Orbit')
state = KeplarianState(orbit, name='Test State')

import logging

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

if __name__ == '__main__':
    orbit.e = 0.8
    state.r = 18500 * ureg.km
    state.ta = 45 * ureg.deg
    print(orbit)
    print()
    print(state)

    plt.figure(1)
    plotting.plot_orbit_fixed(orbit)
    plt.show(block=False)

    pass
