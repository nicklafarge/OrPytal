########### Standard ###########
import itertools
import logging
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames, bodies
from orpytal.common import units

########### External ###########
import numpy as np
from orpytal.planet_constants import BODIES_532, BODIES

logging.basicConfig()
logging.getLogger().setLevel(logging.ERROR)

earth = bodies.earth

possible_values = ['a', 'e', 'rp', 'ra', 'e', 'p', 'h', 'period', 'se']
two_value_pairs = [i for i in itertools.combinations(possible_values, 2)]
print(len(two_value_pairs))

impossible_pairs = [
    ('p', 'h'),         # Direct dependency (h = sqrt(p * mu))
    ('a', 'se'),        # Direct dependency (se = -mu/sqrt(a))
    ('a', 'period'),    # Direct dependency (period = 2pi sqrt(a^3/mu))
    ('period', 'se'),   # Indirect dependency (se <-> a <-> period)
]


class TestOrbitCreation(unittest.TestCase):

    def test_two_value_pairs_elliptic(self):
        orbit = Orbit(earth, a=51000 * units.km, e=0.7)
        for pair in two_value_pairs:
            if pair[0] == pair[1]:
                continue

            if pair in impossible_pairs:
                continue

            test_orbit = Orbit(earth)
            setattr(test_orbit, pair[0], getattr(orbit, pair[0]))
            setattr(test_orbit, pair[1], getattr(orbit, pair[1]))
            same = orbit.compare(test_orbit)
            assert same


if __name__ == '__main__':
    unittest.main()
