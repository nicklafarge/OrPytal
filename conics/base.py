########### Standard ###########

########### Local ###########
from common import ureg, Q_
import frames

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
        if isinstance(value, ureg.Quantity):
            self._value = value.to(self.units)
        elif isinstance(value, frames.Vector):
            if isinstance(value.value, ureg.Quantity):
                self._value = frames.Vector(value.value.to(self.units), value.frame)
            else:
                self._value = frames.Vector(value.value * self.units, value.frame)
        else:
            self._value = value * self.units

        self.evaluated = True

    def check_satisfied(self, obj, req):
        return hasattr(obj, '_' + req) and getattr(obj, '_' + req).evaluated

    def state_orbit_satisfied(self, state, orbit, requirements):
        return all([self.check_satisfied(state, req) or self.check_satisfied(orbit, req) for req in requirements])
