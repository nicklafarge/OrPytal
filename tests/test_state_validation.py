########### Standard ###########
import logging
import itertools
import unittest

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, integration, frames, Trajectory
from orpytal.errors import InvalidInputError

########### External ###########
import numpy as np

orbit = Orbit(bodies.earth,
              a=51000 * units.km,
              e=0.7,
              raan=10 * units.deg,
              arg_periapsis=10 * units.deg,
              inclination=45 * units.deg)


class TestStateValidationMethods(unittest.TestCase):
    def test_validate_position(self):
        st = orbit.get_state()

        st.r = 10 * units.km
        assert st.r == None

        st.r = 1000000 * units.km
        assert st.r == None

        st.r = 61000
        assert st.r == 61000 * units.km

    def test_validate_velocity(self):
        st = orbit.get_state()

        st.v = 1 * units("km/s")
        assert st.v == None

        st.v = 10 * units("km/s")
        assert st.v == None

        st.v = 6
        assert st.v == 6 * units("km/s")


if __name__ == '__main__':
    unittest.main()
