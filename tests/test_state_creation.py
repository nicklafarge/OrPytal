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

possible_eccentric_values = [
    ('ta'),
    ('ttp'),
    ('E'),
    ('arg_latitude'),
    ('r', 'ascending'),
    ('position', 'velocity'),
    ('r', 'fpa')
]

impossible_vals = [
]

unsupported_vals = [

]

a = 51000 * units.km
i = 1.85 * units.deg
raan = 49.562 * units.deg
argp = 286.537 * units.deg
ta = 45 * units.deg

circular_orbit = Orbit(earth, a=a, e=0, raan=raan, arg_periapsis=argp, i=i)
elliptic_orbit = Orbit(earth, a=a, e=0.7, raan=raan, arg_periapsis=argp, i=i)

ascending_elliptic_state = elliptic_orbit.get_state(ta=45 * units.deg)


class TestOrbitCreation(unittest.TestCase):
    def test_ascending_state_matches_poliastro(self):
        elliptic_poliastro = PoliastroOrbit.from_classical(Earth, a.m * u.km, 0.7 * u.one, i.m * u.deg,
                                                           raan.m * u.deg, argp.m * u.deg, ta.m * u.deg)
        ascending_elliptic_state.compare_poliastro(elliptic_poliastro)


if __name__ == '__main__':
    unittest.main()
