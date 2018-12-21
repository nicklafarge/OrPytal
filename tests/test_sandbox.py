########### Standard ###########
import logging
import itertools

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting

########### External ###########
import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

earth = bodies.earth

possible_values = ['a', 'e', 'rp', 'ra', 'e', 'p', 'h', 'period', 'se']
two_value_pairs = [i for i in itertools.combinations(possible_values, 2)]

impossible_pairs = [
    ('p', 'h'),         # Direct dependency (h = sqrt(p * mu))
    ('a', 'se'),        # Direct dependency (se = -mu/sqrt(a))
    ('a', 'period'),    # Direct dependency (period = 2pi sqrt(a^3/mu))
    ('period', 'se'),   # Indirect dependency (se <-> a <-> period)
]

orbit = Orbit(earth, a=51000 * units.km, e=1e-4)
st = KeplarianState(orbit, r=51000*units.km)
