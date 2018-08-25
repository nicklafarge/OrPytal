########### Standard ###########
import unittest

########### Local ###########
from src.orpytal import Orbit, KeplarianState
from common import units

########### External ###########
import numpy as np


class TestOrbit(unittest.TestCase):
    def test_ta(self):
        orbit = Orbit('earth')
        state = KeplarianState(orbit)

        orbit.e = 0.1
        orbit.p = 18400000 * units.m

        state.r = 18500 * units.km

        ta_val = 1.6248767384090719
        assert (np.isclose(state.ta, ta_val))

    def test_from_state(self):
        orbit = Orbit('earth')
        state = KeplarianState(orbit)

        orbit.e = 0.1

        state.r = 18500 * units.km
        state.ta = 90 * units.deg

        orbit.from_state(state)

        assert (np.isclose(orbit.p, 18500.0))
        assert (np.isclose(orbit.h, 85872.6269))


if __name__ == '__main__':
    unittest.main()
