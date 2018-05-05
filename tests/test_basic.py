########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, State
from common import ureg, Q_

########### External ###########
import numpy as np

orbit = Orbit('earth')
state = State(orbit)


class TestOrbit(unittest.TestCase):
    def test_ta(self):
        orbit.e = 0.1
        orbit.p = 18400000 * ureg.m

        state.r = 18500 * ureg.km

        ta_val = 1.6248767384090719
        assert(np.isclose(state.ta, ta_val))

if __name__ == '__main__':
    unittest.main()