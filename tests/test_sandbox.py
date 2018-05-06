########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState
from common import ureg, Q_

########### External ###########
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit('earth', name='Test Orbit')
state = KeplarianState(orbit, name='Test State')

if __name__ == '__main__':
    orbit.e = 0.1

    state.r = 18500 * ureg.km
    state.ta = 90 * ureg.deg

    orbit.from_state(state)

    state.position.orbit_fixed()

    print(orbit)
    print()
    print(state)

    pass
