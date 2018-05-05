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
        self._arg_latitude = ArgumentOfLatitude()

        self.vars = [
            self._r,
            self._ta,
            self._arg_latitude
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
    def arg_latitude(self):
        return self._arg_latitude.value

    @arg_latitude.setter
    def arg_latitude(self, arg_latitude=None):
        self._arg_latitude.value = arg_latitude
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


class ArgumentOfLatitude(StateValue):
    symbol = 'arg_latitude'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [
            ('arg_periapsis', 'ta')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # omega + ta
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.arg_periapsis + state.ta

        # Requirements not met
        else:
            return False

        return True


class RotationMatrixEr(StateValue):
    symbol = 'dcm_er'

    def __init__(self):
        super().__init__(ureg.dimensionless)
        self.orbit_requirements = [
            TrueAnomaly.symbol
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = np.array(
                [[np.cos(state.ta), np.sin(state.ta), 0],
                 [-np.sin(state.ta), np.cos(state.ta), 0],
                 [0, 0, 1]]
            )

        # Requirements not met
        else:
            return False

        return True


class RotationMatrixRv(StateValue):
    symbol = 'dcm_rv'

    def __init__(self):
        super().__init__(ureg.dimensionless)
        self.orbit_requirements = [
            'fpa'
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = np.array(
                [[np.sin(state.fpa), 0, np.cos(state.fpa)],
                 [np.cos(state.fpa), 0, -np.sin(state.fpa)],
                 [0, 1, 0]])

        # Requirements not met
        else:
            return False

        return True


class RotationMatrixRi(StateValue):
    symbol = 'dcm_ri'

    def __init__(self):
        super().__init__(ureg.dimensionless)
        self.orbit_requirements = [
            'arg_latitude', 'ascending_node', 'inclination'
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            Omega = orbit.ascending_node
            theta = state.arg_latitude
            i = orbit.inclination
            dcm_ri = np.zeros((3, 3))
            dcm_ri[0, 0] = np.cos(Omega) * np.cos(theta) - np.sin(Omega) * np.cos(i) * np.sin(theta)
            dcm_ri[0, 1] = -np.cos(Omega) * np.sin(theta) - np.sin(Omega) * np.cos(i) * np.cos(theta)
            dcm_ri[0, 2] = np.sin(Omega) * np.sin(i)
            dcm_ri[1, 0] = np.sin(Omega) * np.cos(theta) + np.cos(Omega) * np.cos(i) * np.sin(theta)
            dcm_ri[1, 1] = -np.sin(Omega) * np.sin(theta) + np.cos(Omega) * np.cos(i) * np.cos(theta)
            dcm_ri[1, 2] = -np.cos(Omega) * np.sin(i)
            dcm_ri[2, 0] = np.sin(i) * np.sin(theta)
            dcm_ri[2, 1] = np.sin(i) * np.cos(theta)
            dcm_ri[2, 2] = np.cos(i)

            self.value = dcm_ri

        # Requirements not met
        else:
            return False

        return True
