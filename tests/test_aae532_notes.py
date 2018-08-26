########### Standard ###########
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState,frames
from orpytal.common import units
from orpytal.planet_constants import BODIES_532

########### External ###########
import numpy as np

earth = BODIES_532['EARTH']
earth.radius = 6378.137 * units.km

class TestOrbit(unittest.TestCase):

    def test_JS_3Dex1(self):
        orbit = Orbit(earth, name='Test Orbit')
        state = KeplarianState(orbit, name='Test State')

        orbit.a = 8 * earth.radius
        orbit.e = 0.7
        orbit.i = 30 * units.deg
        orbit.ascending_node = 60 * units.deg
        orbit.arg_periapsis = 90 * units.deg
        state.ta = 90 * units.deg

        state.position.inertial()
        pass

    def test_Iex(self):
        orbit = Orbit(earth, name='Test Orbit')
        state = KeplarianState(orbit, name='Test State')

        state.position = frames.Vector(orbit, state, np.array([1.6772, -1.6722, 2.3719]) * earth.radius, frames.InertialFrame)
        state.velocity = frames.Vector(orbit, state, np.array([3.1574, 2.4987, 0.4658]), frames.InertialFrame)


        pass


if __name__ == '__main__':
    unittest.main()
