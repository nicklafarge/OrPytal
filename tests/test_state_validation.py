########### Standard ###########
import logging
import itertools
import unittest

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, integration, frames, Trajectory
from orpytal.errors import InvalidInputError

########### External ###########
import numpy as np

logging.disable(logging.CRITICAL)

orbit = Orbit(bodies.earth,
              a=51000 * units.km,
              e=0.7,
              raan=10 * units.deg,
              arg_periapsis=10 * units.deg,
              inclination=45 * units.deg)


class TestStateValidationMethods(unittest.TestCase):
    """
        Test state validation methods (that don't allow you to set wrong balues
    """

    def test_validate_position(self):
        """
            Test position validation (can't be inside periapsis, etc)
        """
        st = orbit.get_state()

        st.r = 10 * units.km
        assert st.r == None

        st.r = 1000000 * units.km
        assert st.r == None

        st.r = 61000
        assert st.r == 61000 * units.km

    def test_validate_velocity(self):
        """
            Test velocity validation (can't be faster than v at rp)
        """
        st = orbit.get_state()

        st.v = 1 * units("km/s")
        assert st.v == None

        st.v = 10 * units("km/s")
        assert st.v == None

        st.v = 6
        assert st.v == 6 * units("km/s")

    def test_constructor_syntax(self):
        """
            Test validation from KeplarianState constructor
        """
        st = KeplarianState(orbit, r=10 * units.km)
        assert st.r == None

    def test_get_state_syntax(self):
        """
            Test validation from get_state syntax
        """
        st = orbit.get_state(r=10 * units.km)
        assert st.r == None

    def test_prohibit_resetting_variable(self):
        """
            Test that if a variable is already set, you can't set it again
        """
        bad_r = 14000*units.km
        st = orbit.get_state(ta=60)
        st.r = bad_r
        assert st.r != bad_r


if __name__ == '__main__':
    unittest.main()
