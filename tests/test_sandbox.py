########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState, plotting
from common import units, Q_

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
    orbit.p = 37800 * units.km
    state.ta = -175 * units.deg
    print(orbit)
    print()
    print(state)

    plt.figure(1)
    plotting.plot_orbit(orbit)
    plotting.plot_state(state)
    plt.show(block=False)

    plt.figure(2)
    plotting.animate_orbit_fixed(orbit)

    pass
