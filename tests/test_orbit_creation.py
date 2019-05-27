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

possible_values = ['a', 'e', 'rp', 'ra', 'e', 'p', 'h', 'period', 'se', 'b']
two_value_pairs = [i for i in itertools.combinations(possible_values, 2)]
print(len(two_value_pairs))

impossible_pairs = [
    ('p', 'h'),  # Direct dependency (h = sqrt(p * mu))
    ('a', 'se'),  # Direct dependency (se = -mu/sqrt(a))
    ('a', 'period'),  # Direct dependency (period = 2pi sqrt(a^3/mu))
    ('period', 'se'),  # Indirect dependency (se <-> a <-> period)
]

# I don't know if these are impossible
unsupported_pairs = [
    ('rp', 'period'),  # I suspect ths isn't possible
    ('p', 'period'),  # I suspect ths isn't possible
    ('h', 'period'),  # I suspect ths isn't possible
    ('e', 'b'),  # Should be possible
    ('rp', 'b'),  # Should be possible
    ('ra', 'b'),  # Should be possible
    ('p', 'b'),  # Should be possible
    ('h', 'b'),  # Should be possible
    ('period', 'b'),  # Should be possible
]

a = 51000 * units.km
i = 1.85 * units.deg
raan = 49.562 * units.deg
argp = 286.537 * units.deg
ta = 23.33 * units.deg

circular_orbit = Orbit(earth, a=a, e=0, raan=raan, arg_periapsis=argp, i=i)
elliptic_orbit = Orbit(earth, a=a, e=0.7, raan=raan, arg_periapsis=argp, i=i)

circular_state = circular_orbit.get_state(ta=0)
elliptic_state = elliptic_orbit.get_state(ta=0)

to = Orbit(earth, raan=raan, arg_periapsis=argp, i=i)
to.a = a
to.e = 0.2
x = 1


class TestOrbitCreation(unittest.TestCase):
    def test_angle_units_after_orbit_creation_regular(self):
        """
        Check that angles have proper units (and haven't been converted to dimensionless) for regular orbit creation
        """
        orbit_test = Orbit(earth)
        orbit_test.i = i
        assert orbit_test.i.units == units.rad

    def test_angle_units_after_orbit_creation_constructor(self):
        """
        Check that angles have proper units (and haven't been converted to dimensionless) for constructor orbit creation
        """
        orbit_test = Orbit(earth, i=i)
        assert orbit_test.i.units == units.rad

    def test_circular_orbit_matches_poliastro(self):
        """
        Check that the circular orbit here is the same as if created using poliastro
        """

        circular_poliastro = PoliastroOrbit.from_classical(Earth, a.m * u.km, 0 * u.one, i.m * u.deg,
                                                           raan.m * u.deg, argp.m * u.deg, ta.m * u.deg)

        circular_orbit.compare_poliastro(circular_poliastro)

    def test_elliptic_orbit_matches_poliastro(self):
        """
        Check that the elliptic orbit here is the same as if created using poliastro
        """

        elliptic_poliastro = PoliastroOrbit.from_classical(Earth, a.m * u.km, 0.7 * u.one, i.m * u.deg,
                                                           raan.m * u.deg, argp.m * u.deg, ta.m * u.deg)
        elliptic_orbit.compare_poliastro(elliptic_poliastro)

    def test_two_value_pairs_elliptic(self):
        """
        Check that all supported two value pairs yield the same orbit as all other two value pairs (elliptic orbit)
        """
        for pair in two_value_pairs:
            if pair[0] == pair[1]:
                continue

            # Skip over unsupported pairs
            if pair in impossible_pairs or pair in unsupported_pairs:
                continue

            test_orbit = Orbit(earth, i=i, raan=raan, arg_periapsis=argp)
            setattr(test_orbit, pair[0], getattr(elliptic_orbit, pair[0]))
            setattr(test_orbit, pair[1], getattr(elliptic_orbit, pair[1]))
            same = elliptic_orbit.compare(test_orbit)
            if not same:
                logging.error(pair)
            assert same
            assert test_orbit.type() == OrbitType.Elliptic
            assert not test_orbit.circular()

    def test_two_value_pairs_circular(self):
        """
        Check that all supported two value pairs yield the same orbit as all other two value pairs (circular orbit)
        """

        for pair in two_value_pairs:
            if pair[0] == pair[1]:
                continue

            # Skip over unsupported pairs
            if pair in impossible_pairs or pair in unsupported_pairs:
                continue

            test_orbit = Orbit(earth, i=i, raan=raan, arg_periapsis=argp)
            setattr(test_orbit, pair[0], getattr(circular_orbit, pair[0]))
            setattr(test_orbit, pair[1], getattr(circular_orbit, pair[1]))
            same = circular_orbit.compare(test_orbit)
            if not same:
                logging.error(pair)
            assert same
            assert test_orbit.type() == OrbitType.Elliptic
            assert test_orbit.circular()

    def test_not_ascending(self):
        orbit = Orbit(earth, a=a, e=0.4)
        st = orbit.get_state(r=49000)
        assert st._ascending == None
        assert st.ta == None


if __name__ == '__main__':
    unittest.main()
