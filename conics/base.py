########### Standard ###########

########### Local ###########
from common import units, Q_
import frames
import conics_utils

########### External ###########
import numpy as np
import scipy as sp


class OrbitBase(object):
    symbol = ''

    def __init__(self, units):
        self.units = units

        self.requirements = None
        self._value = None
        self.evaluated = False

    def __str__(self):
        return '{}: {}'.format(self.symbol, self._value)

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value=None):
        if isinstance(value, units.Quantity):
            self._value = value.to(self.units)
        elif isinstance(value, frames.Vector):
            if isinstance(value.value, units.Quantity):
                value.value = value.value.to(self.units)
                self._value = value
            else:
                value.value = value.value * self.units
                self._value = value
        else:
            self._value = value * self.units

        # Resolve signs for units
        if self.units == units.rad:
            self._value = self._value % (2* np.pi)

        self.evaluated = True

    def check_satisfied(self, obj, req):
        return conics_utils.check_satisfied(obj, req)

    def state_orbit_satisfied(self, state, orbit, requirements):
        return conics_utils.state_orbit_satisfied(state, orbit, requirements)