########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState
from common import ureg, Q_

########### External ###########
import numpy as np


class TestOrbit(unittest.TestCase):
    def test_ta(self):
        orbit = Orbit('earth')
        state = KeplarianState(orbit)

        orbit.e = 0.1
        orbit.p = 18400000 * ureg.m

        state.r = 18500 * ureg.km

        ta_val = 1.6248767384090719
        assert (np.isclose(state.ta, ta_val))

    def test_from_state(self):
        orbit = Orbit('earth')
        state = KeplarianState(orbit)

        orbit.e = 0.1

        state.r = 18500 * ureg.km
        state.ta = 90 * ureg.deg

        orbit.from_state(state)

        assert (np.isclose(orbit.p, 18500.0))
        assert (np.isclose(orbit.h, 85872.6269))


if __name__ == '__main__':
    unittest.main()
