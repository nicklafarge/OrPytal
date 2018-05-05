########### Standard ###########

########### Local ###########
from common import ureg, Q_

########### External ###########
import numpy as np
import scipy as sp


class OrbitBase(object):
    def __init__(self, units):
        self.units = units

        self.requirements = None
        self._value = None
        self.evaluated = False

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value=None):
        if isinstance(value, ureg.Quantity):
            self._value = value.to(self.units)
        else:
            self._value = value * self.units

        self.evaluated = True

    def check_satisfied(self, obj, req):
        return hasattr(obj, '_' + req) and getattr(obj, '_' + req).evaluated

    def state_orbit_satisfied(self, state, orbit, requirements):
        return all([self.check_satisfied(state, req) or self.check_satisfied(orbit, req) for req in requirements])
