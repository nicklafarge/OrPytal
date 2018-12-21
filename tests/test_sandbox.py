########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting

########### External ###########
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

earth = bodies.earth

possible_values = ['a', 'e', 'rp', 'ra', 'e', 'p', 'h', 'period', 'se']
two_value_pairs = [i for i in itertools.combinations(possible_values, 2)]

impossible_pairs = [
    ('p', 'h'),         # Direct dependency (h = sqrt(p * mu))
    ('a', 'se'),        # Direct dependency (se = -mu/sqrt(a))
    ('a', 'period'),    # Direct dependency (period = 2pi sqrt(a^3/mu))
    ('period', 'se'),   # Indirect dependency (se <-> a <-> period)
]

orbit = Orbit(earth, a=51000 * units.km, e=0)

for pair in two_value_pairs:
    if pair[0] == pair[1]:
        continue

    if pair in impossible_pairs:
        continue

    test_orbit = Orbit(earth)
    setattr(test_orbit, pair[0], getattr(orbit, pair[0]))
    setattr(test_orbit, pair[1], getattr(orbit, pair[1]))
    same = orbit.compare(test_orbit)
    if not same:
        print(pair)
