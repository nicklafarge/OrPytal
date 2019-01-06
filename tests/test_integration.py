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
    def test_numerical_intergration(self):
        traj = orbit.propagate_orbit()
        assert np.isclose(orbit.rp, traj.end().r)


if __name__ == '__main__':
    unittest.main()
