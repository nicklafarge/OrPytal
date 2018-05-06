########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, KeplarianState
from common import ureg, Q_
import frames

########### External ###########
import numpy as np
from planet_constants import BODIES_532

earth = BODIES_532['EARTH']
earth.radius = 6371 * ureg.km


class TestOrbit(unittest.TestCase):

    def test_ps3_problem_3(self):
        orbit = Orbit(earth, name='Test Orbit')
        state = KeplarianState(orbit, name='Test State')

        state.r = orbit.central_body.radius + 2000 * ureg.km
        vel_vector = np.array([-1.2, 6.7, 0.0]) * ureg.km / ureg.s
        state.velocity = frames.Vector(orbit, state, vel_vector, frames.RotatingFrame)

        print(orbit)
        print()
        print(state)

        # Assert Orbit Quantities
        assert (np.isclose(orbit.p, 7891.62633))
        assert (np.isclose(orbit.e, 0.17829470376))
        assert (np.isclose(orbit.h, 56085.7))
        assert (np.isclose(orbit.a, 8150.7299))
        assert (np.isclose(orbit.b, 7648.4885248))
        assert (np.isclose(orbit.period, 7323.281580402567))
        assert (np.isclose(orbit.n, 0.0008579740159))
        assert (np.isclose(orbit.se, -24.451824967))
        assert (np.isclose(orbit.rp, 9603.961874))
        assert (np.isclose(orbit.ra, 6697.497928))

        # Assert State Quanitites
        assert (np.isclose(state.ta, -1.8977792920))
        assert (np.isclose(state.v, 6.80661443))
        assert (np.isclose(state.fpa, -0.17722538))
        assert (np.isclose(state.t_since_rp, -1802.7593))
        assert (np.isclose(state.M, -1.546720643))
        assert (np.isclose(state.E, -1.7229553535))

        pass


if __name__ == '__main__':
    unittest.main()
