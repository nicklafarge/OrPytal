########### Standard ###########
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames
from orpytal.common import units

########### External ###########
import numpy as np
from orpytal.planet_constants import BODIES_532, BODIES

earth = BODIES_532['EARTH']


class TestOrbit(unittest.TestCase):

    def test_ps3_problem_3(self):
        earth.radius = 6378.137 * units.km
        orbit = Orbit(earth, name='Test Orbit')
        state = KeplarianState(orbit, name='Test State')

        state.r = orbit.central_body.radius + 2000 * units.km
        vel_vector = np.array([-1.2, 6.7, 0.0]) * units.km / units.s
        state.velocity = vel_vector, frames.RotatingFrame

        print(orbit)
        print()
        print(state)

        from orpytal import plotting
        from matplotlib import pyplot as plt
        plotting.plot_orbit(orbit, frames.OrbitFixedFrame)
        plotting.plot_state(state, frames.OrbitFixedFrame)
        plt.show(block=False)

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

        # Assert State Quantities
        assert (np.isclose(state.ta, -1.8977792920))
        assert (np.isclose(state.v, 6.80661443))
        assert (np.isclose(state.fpa, -0.17722538))
        assert (np.isclose(state.t_since_rp, -1802.7593))
        assert (np.isclose(state.M, -1.546720643))
        assert (np.isclose(state.E, -1.7229553535))

        pass

    def test_ps4_problem_4(self):
        earth.radius = 6378.137 * units.km
        orbit = Orbit(earth, name='Test Orbit')
        state = KeplarianState(orbit, name='Test State')
        orbit.a = 10 * earth.radius
        orbit.e = 0.4
        orbit.i = 30 * units.deg
        orbit.ascending_node = 45 * units.deg
        orbit.arg_periapsis = -90 * units.deg
        state.ta = 135 * units.deg

        assert (np.isclose(state.r, 74707))
        assert (np.isclose(state.v, 2.10276))
        assert (np.isclose(orbit.h, 146135))
        assert (np.isclose(state.fpa.to('deg'), 21.52399*units.deg))
        assert (np.isclose(state.E.to('deg'), 115.355*units.deg))


        pass


if __name__ == '__main__':
    unittest.main()
