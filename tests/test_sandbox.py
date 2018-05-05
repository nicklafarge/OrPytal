########### Standard ###########
import unittest

########### Local ###########
from conics import Orbit, State
from common import ureg, Q_

########### External ###########
import numpy as np

orbit = Orbit('earth')
state = State(orbit)

if __name__ == '__main__':
    orbit.e = 0.1

    state.r = 18500 * ureg.km
    state.ta = 90 * ureg.deg

    orbit.from_state(state)

