########### Standard ###########
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames, bodies, OrbitType
from orpytal.common import units
from orpytal.planet_constants import CentralBody

########### External ###########
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit as PoliastroOrbit
from astropy import units as u

# Create the central body
earth_poliastro = CentralBody(
    name='Earth (Poliastro)',
    radius=bodies.earth.radius,
    mu=Earth.k.value / 1e9
)

earth = earth_poliastro

# Create test orbit/state objects for the tests
hyperbolic_orbit = Orbit(earth, a=-51000, e=1.1, i=0.1, raan=0.4, arg_periapsis=1.1)
ascending_state = KeplarianState(hyperbolic_orbit, ta=30 * units.deg)


class TestHyperbolicOrbitCreation(unittest.TestCase):
    """
        Test creation and methods pertaining to hyperbolic orbits
    """

    def test_hyperbolic_orbit_matches_poliastro(self):
        """
        Check that the hyperbolic orbit here is the same as if created using poliastro
        """

        hyperoblic_poliastro = PoliastroOrbit.from_classical(Earth,
                                                             hyperbolic_orbit.a.m * u.km,
                                                             hyperbolic_orbit.e.m * u.one,
                                                             hyperbolic_orbit.i.m * u.rad,
                                                             hyperbolic_orbit.raan.m * u.rad,
                                                             hyperbolic_orbit.arg_periapsis.m * u.rad,
                                                             ascending_state.ta.m * u.rad)

        hyperbolic_orbit.compare_poliastro(hyperoblic_poliastro)

    def test_orbit_type(self):
        """
            Check that the orbit here is correctly identified as hyperbolic
        """

        assert (hyperbolic_orbit.type() == OrbitType.Hyperbolic)

    def test_creation_from_position_and_velocity_hyperbolic(self):
        """
        Test creation from inertial position and velocity for a hyperbolic orbit
        """

        test_orbit = Orbit(earth,
                           i=hyperbolic_orbit.i,
                           raan=hyperbolic_orbit.raan,
                           arg_periapsis=hyperbolic_orbit.arg_periapsis)
        new_st = KeplarianState(test_orbit)
        new_st.position = ascending_state.position
        new_st.velocity = ascending_state.velocity

        assert hyperbolic_orbit.compare(test_orbit)
        assert test_orbit.type() == OrbitType.Hyperbolic


if __name__ == '__main__':
    unittest.main()
