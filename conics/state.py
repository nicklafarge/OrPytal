########### Standard ###########

########### Local ###########
from base import OrbitBase
from common import ureg, Q_
import frames

########### External ###########
import numpy as np
import scipy as sp


class State(object):
    def __init__(self, orbit, name=''):
        self.orbit = orbit
        self.name = name

        self._ascending = True

        self._r = PositionMagnitude()
        self._ta = TrueAnomaly()
        self._arg_latitude = ArgumentOfLatitude()
        self._v = VelocityMagnitude()
        self._pos = PositionVector()
        self._vel = VelocityVector()
        self._fpa = FlightPathAngle()
        self._t_since_rp = TimeSincePeriapsis()
        self._M = MeanAnomaly()
        self._E = EccentricAnomaly()

        self.vars = [
            self._r,
            self._ta,
            self._arg_latitude,
            self._v,
            self._pos,
            self._vel,
            self._fpa,
            self._t_since_rp,
            self._M,
            self._E
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
    def v(self):
        return self._v.value

    @v.setter
    def v(self, v=None):
        self._v.value = v
        self.set_vars()

    @property
    def pos(self):
        return self._pos.value

    @pos.setter
    def pos(self, pos=None):
        self._pos.value = pos
        self.set_vars()

    @property
    def vel(self):
        return self._vel.value

    @vel.setter
    def vel(self, vel=None):
        self._vel.value = vel
        self.set_vars()

    @property
    def fpa(self):
        return self._fpa.value

    @fpa.setter
    def fpa(self, fpa=None):
        self._fpa.value = fpa
        self.set_vars()

    @property
    def M(self):
        return self._M.value

    @M.setter
    def M(self, M=None):
        self._M.value = M
        self.set_vars()

    @property
    def E(self):
        return self._E.value

    @E.setter
    def E(self, E=None):
        self._E.value = E
        self.set_vars()

    @property
    def t_since_rp(self):
        return self._t_since_rp.value

    @t_since_rp.setter
    def t_since_rp(self, t_since_rp=None):
        self._t_since_rp.value = t_since_rp
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


class PositionMagnitude(StateValue):
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


class FlightPathAngle(StateValue):
    symbol = 'fpa'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [
            ('h', 'r', 'v'),
            ('velocity')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # acos(h/(rv))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = state.ascending_sign() * np.arccos(orbit.h / (state.r * state.v))

        # atan(vr/vt)
        if self.satisfied(state, orbit, self.orbit_requirements[1]):
            fpa = np.arctan(state.velocity[0] / state.velocity[1])
            self.value = state.angle_check_tan(fpa)

        # Requirements not met
        else:
            return False

        return True


class VelocityMagnitude(StateValue):
    symbol = 'v'

    def __init__(self):
        super().__init__(ureg.radian)
        self.orbit_requirements = [
            ('r', 'se'),
            ('r', 'a'),
            ('velocity')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # sqrt( 2(se + mu/r))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(2 * (orbit.se + orbit.central_body.mu / state.r))
        # sqrt( mu(2/r-1/a) )
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = np.sqrt(orbit.central_body.mu * (2. / state.r - 1 / orbit.a))
        # |v|
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            self.value = state.velocity.norm()

        # Requirements not met
        else:
            return False

        return True


class PositionVector(StateValue):
    symbol = 'position'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('r')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # r (rhat)
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = frames.Vector(state.r, frames.RotatingFrame)

        # Requirements not met
        else:
            return False

        return True


class VelocityVector(StateValue):
    symbol = 'velocity'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('ta', 'e', 'h', 'ta'),
            ('fpa' 'v')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # Rotating Frame
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            vr = (orbit.central_body.mu * orbit.e) / orbit.h * np.sin(state.ta)
            vt = orbit.central_body.mu / orbit.h * (1 + orbit.e * np.cos(state.ta))
            self.value = frames.Vector(np.array([vr, vt, 0]), frames.RotatingFrame)

        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = frames.Vector(np.array([state.v * np.sin(state.fpa), state.v * np.cos(state.fpa), 0]),
                                       frames.RotatingFrame)
        # Requirements not met
        else:
            return False

        return True


class TimeSincePeriapsis(StateValue):
    symbol = 't_since_rp'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('E', 'e', 'n'),
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # Rotating Frame
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = (state.E - orbit.e * np.sin(state.E)) / orbit.n

        # Requirements not met
        else:
            return False

        return True


class MeanAnomaly(StateValue):
    symbol = 'M'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('n', 't_since_rp'),
            ('E', 'e')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # n * ttp
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.n * state.t_since_rp

        # E - e sin(E)
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = state.E - orbit.e * np.sin(state.E)

        # Requirements not met
        else:
            return False

        return True


class EccentricAnomaly(StateValue):
    symbol = 'E'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('a', 'r', 'e'),
            ('e', 'M')
        ]

    def set(self, state, orbit):
        if self.evaluated:
            return False

        # acos((a-r)/(ae))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = state.ascending_sign() * np.arccos(
                (orbit.a - state.r) / (orbit.a * orbit.e))
        if self.satisfied(state, orbit, self.orbit_requirements[1]):
            self._iterative_eccentric_anomaly(state)

        # Requirements not met
        else:
            return False

        return True

    def _iterative_eccentric_anomaly(self, state, orbit, **kwargs):
        return self._find_eccentric_anomaly_recursively(orbit, state.M, state.M, **kwargs)

    def _find_eccentric_anomaly_recursively(self, orbit, E, M, tol=1e-12):

        E1 = E - (E - orbit.e * np.sin(E) - M) / (1 - orbit.e * np.cos(E))
        dE = np.fabs(E - E1)

        if dE < tol:
            return E

        return self._find_eccentric_anomaly_recursively(orbit, E1, M, tol=tol)
