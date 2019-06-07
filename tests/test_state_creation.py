########### Standard ###########
import itertools
import logging
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames, bodies, OrbitType
from orpytal.common import units
from orpytal.planet_constants import CentralBody

########### External ###########
import numpy as np
import poliastro
from poliastro.bodies import Earth
from poliastro.twobody import Orbit as PoliastroOrbit
from astropy import units as u

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

earth_poliastro = CentralBody(
    name='Earth (Poliastro)',
    radius=bodies.earth.radius,
    mu=Earth.k.value / 1e9
)

earth = earth_poliastro

# Create some test values
a = 51000 * units.km
i = 1.85 * units.deg
raan = 49.562 * units.deg
argp = 286.537 * units.deg
ta = 45 * units.deg

# Create test OrPytal orbits
circular_orbit = Orbit(earth, a=a, e=0, raan=raan, arg_periapsis=argp, i=i)
elliptic_orbit = Orbit(earth, a=a, e=0.7, raan=raan, arg_periapsis=argp, i=i)

# Create test OrPytal states
ascending_elliptic_state = elliptic_orbit.get_state(ta=45 * units.deg)
descending_elliptic_state = elliptic_orbit.get_state(ta=-45 * units.deg)


class TestOrbitCreation(unittest.TestCase):
    """
        Test Keplarian state creation
    """
    def test_ascending_state_matches_poliastro(self):
        """
            Test that the sample ascending state parameters match those in poliastro
        """
        elliptic_poliastro = PoliastroOrbit.from_classical(Earth, a.m * u.km, 0.7 * u.one, i.m * u.deg,
                                                           raan.m * u.deg, argp.m * u.deg, ta.m * u.deg)
        ascending_elliptic_state.compare_poliastro(elliptic_poliastro)

    def test_descending_state_matches_poliastro(self):
        """
            Test that the sample descending state parameters match those in poliastro
        """
        elliptic_poliastro = PoliastroOrbit.from_classical(Earth, a.m * u.km, 0.7 * u.one, i.m * u.deg,
                                                           raan.m * u.deg, argp.m * u.deg,
                                                           descending_elliptic_state.ta.m * u.rad)
        descending_elliptic_state.compare_poliastro(elliptic_poliastro)

    def test_ascending_state_creation(self):
        """
            Test that various ascending creation state methods yield the same result
        """
        true_anomaly_state = ascending_elliptic_state

        # Bigger TA
        assert true_anomaly_state.compare(elliptic_orbit.get_state(ta=true_anomaly_state.ta + 2 * np.pi))

        # Time since Periapsis
        assert true_anomaly_state.compare(elliptic_orbit.get_state(t_since_rp=true_anomaly_state.t_since_rp))

        # Eccentric Anomaly
        assert true_anomaly_state.compare(elliptic_orbit.get_state(E=true_anomaly_state.E))

        # ascending r
        assert true_anomaly_state.compare(elliptic_orbit.get_state(ascending=True, r=true_anomaly_state.r))

        # r,fpa
        assert true_anomaly_state.compare(elliptic_orbit.get_state(fpa=true_anomaly_state.fpa, r=true_anomaly_state.r))

    def test_descscending_state_creation(self):
        """
            Test that various descending creation state methods yield the same result
        """
        true_anomaly_state = descending_elliptic_state

        # Bigger TA
        assert true_anomaly_state.compare(elliptic_orbit.get_state(ta=true_anomaly_state.ta - 2 * np.pi))

        # Time since Periapsis
        assert true_anomaly_state.compare(elliptic_orbit.get_state(t_since_rp=true_anomaly_state.t_since_rp))

        # Eccentric Anomaly
        assert true_anomaly_state.compare(elliptic_orbit.get_state(E=true_anomaly_state.E))

        # ascending r
        assert true_anomaly_state.compare(elliptic_orbit.get_state(ascending=False, r=true_anomaly_state.r))

        # r,fpa
        assert true_anomaly_state.compare(elliptic_orbit.get_state(fpa=true_anomaly_state.fpa, r=true_anomaly_state.r))

if __name__ == '__main__':
    unittest.main()
