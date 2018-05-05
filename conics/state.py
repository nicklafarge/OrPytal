########### Standard ###########

########### Local ###########
from base import OrbitBase
from common import ureg, Q_

########### External ###########
import numpy as np
import scipy as sp


class State(object):
    def __init__(self, orbit, name=''):
        self.orbit = orbit
        self.name = name

        self._r = Position()
        self._ta = TrueAnomaly()

        self._ascending = False

        self.vars = [
            self._r, self._ta
        ]

    def set_vars(self):
        for var in self.vars:
            new_value_set = var.set(self, self.orbit)
            if new_value_set:
                self.set_vars()

    @property
    def r(self):
        return self._r.value

    @r.setter
    def r(self, r=None):
        self._r.value = r
        self.set_vars()

    @property
    def ta(self):
        return self._ta.value

    @ta.setter
    def ta(self, ta=None):
        self._ta.value = ta
        self.set_vars()

    @property
    def ascending_sign(self):
        return 1 if self.is_ascending() else -1

    def is_ascending(self):

        if self._ta.evaluated:
            angle_to_check = self.ta
        else:
            return True

        return 0 <= angle_to_check <= np.pi

    def angle_check_tan(self, tan_val):
        if (self.is_ascending() and tan_val > 0) or (not self.is_ascending() and tan_val < 0):
            return tan_val
        else:
            return np.pi + tan_val


class StateValue(OrbitBase):
    def set(self, state, orbit):
        raise NotImplementedError()

    def satisfied(self, state, orbit, requirements):
        return self.state_orbit_satisfied(state, orbit, requirements)


class TrueAnomaly(StateValue):
    symbol = 'ta'

    def __init__(self):
        super().__init__(ureg.radian)
        self.orbit_requirements = [
            ('e', 'p', 'r'),
            ('e', 'E'),
            ('r', 'v', 'fpa')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # acos((1/e) (p/r - 1))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            cos_val = (1. / orbit.e) * (orbit.p / state.r - 1.)
            self.value = state.ascending_sign * np.arccos(cos_val)

        # acos(2 atan(sqrt( (1+e)/(1-e) tan(E/2) )))
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            ta = 2 * np.arctan(np.sqrt((1 + orbit.e) / (1 - orbit.e)) * np.tan(state.E / 2))
            self.value = state.angle_check_tan(ta)

        # atan( ((rv2)/mu cos(g)sin(g)) / ((rv2)/mu cos^2(g)-1) )
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            t = (state.r * state.v ** 2) / orbit.central_body.mu
            ta = np.arctan((t * np.cos(state.fpa) * np.sin(state.fpa)) /
                           (t * np.cos(state.fpa) ** 2 - 1))
            self.value = state.angle_check_tan(ta)

        # Requirements not met
        else:
            return False

        return True


class Position(StateValue):
    symbol = 'r'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('p', 'e', 'ta'),
            ('a', 'e', 'E'),
            ('a', 'E', 'H')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # p/1+ecos(ta)
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1.0 + orbit.e * np.cos(state.ta))

        # a (1-ecosE)
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 - orbit.e * np.cos(state.E))

        # |a|(E cosh H - 1)
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            self.value = abs(orbit.a) * (state.E * np.cosh(state.H) - 1)

        # Requirements not met
        else:
            return False

        return True
