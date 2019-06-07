########### Standard ###########
import logging
import unittest

########### Local ###########
from orpytal import units, Orbit, bodies

import numpy as np

orbit = Orbit(bodies.earth,
              a=51000 * units.km,
              e=0.7,
              raan=10 * units.deg,
              arg_periapsis=10 * units.deg,
              inclination=45 * units.deg)


class TestIntegration(unittest.TestCase):
    """
    Tests for orbit propagation/integration
    """

    def test_numerical_intergration(self):
        """
        Test that, when numerically propagated for one period, the end state is periapsis 
        """
        traj = orbit.propagate_orbit()
        assert np.isclose(orbit.rp, traj.end().r)

    def test_analytic_intergration(self):
        """
        Test that, when analytically propagated for one period, the end state is periapsis 
        """
        traj = orbit.analytic_propagate_full_orbit()
        assert np.isclose(orbit.rp, traj.end().r)


if __name__ == '__main__':
    unittest.main()
