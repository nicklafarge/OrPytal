########### Standard ###########
import logging

########### Local ###########
from orpytal.base import OrbitBase
from orpytal.common import units, conics_utils, orbit_setter, attribute_setter
from orpytal.errors import ParameterUnavailableError
from orpytal import frames

########### External ###########
import numpy as np


class KeplarianState(object):
    def __init__(self, orbit, name=''):
        self.orbit = orbit
        self.name = name

        self._r = PositionMagnitude()
        self._ta = TrueAnomaly()
        self._arg_latitude = ArgumentOfLatitude()
        self._v = VelocityMagnitude()
        self._position = PositionVector()
        self._velocity = VelocityVector()
        self._fpa = FlightPathAngle()
        self._t_since_rp = TimeSincePeriapsis()
        self._M = MeanAnomaly()
        self._E = EccentricAnomaly()

        self._is_ascending = None

        self.vars = [
            self._r,
            self._ta,
            self._arg_latitude,
            self._v,
            self._position,
            self._velocity,
            self._fpa,
            self._t_since_rp,
            self._M,
            self._E
        ]

    def __str__(self):
        x = ['{} Orbit State {}'.format(self.orbit.name, self.name)]
        for var in self.vars:
            if var.evaluated:
                x.append(str(var))

        return '\n'.join(x)

    def set_vars(self):
        for var in self.vars:
            new_value_set = var.set(self, self.orbit)
            if new_value_set:
                logging.debug('Set {} to {}'.format(var.symbol, var.value))
                self.set_vars()
                break

        orbit_changed = self.orbit.from_state(self)
        if orbit_changed:
            self.set_vars()

    def compare(self, other_state):
        # Compare which vars are evaluated
        evaluated_vars = sorted([v.symbol for v in self.vars if v.evaluated])
        other_evaluated_vars = sorted([v.symbol for v in other_state.vars if v.evaluated])
        for i, var in enumerate(evaluated_vars):
            assert var == other_evaluated_vars[i]

        # Compare values of evaluated variables
        for var in self.vars:
            if var.evaluated and hasattr(other_state, var.symbol):
                if isinstance(var.value, units.Quantity):
                    same = np.isclose(var.value, getattr(other_state, var.symbol))
                elif isinstance(var.value, frames.Vector):

                    vec1 = var.value
                    vec2 = getattr(other_state, var.symbol)
                    same = vec1 == vec2
                else:
                    assert False

            if same:
                logging.debug('Checked {} [✓]'.format(var.symbol))
            else:
                logging.warning('Error Found for {} [x]'.format(var.symbol))
    @property
    def r(self):
        return self._r.value

    @r.setter
    @attribute_setter
    def r(self, r):
        self._r.value = r

    @property
    def ta(self):
        return self._ta.value

    @ta.setter
    @attribute_setter
    def ta(self, ta):
        self._ta.value = ta

    @property
    def arg_latitude(self):
        return self._arg_latitude.value

    @arg_latitude.setter
    @attribute_setter
    def arg_latitude(self, arg_latitude):
        self._arg_latitude.value = arg_latitude

    @property
    def v(self):
        return self._v.value

    @v.setter
    @attribute_setter
    def v(self, v):
        self._v.value = v

    @property
    def position(self):
        return self._position.value

    @position.setter
    @attribute_setter
    def position(self, position):
        self._position.value = position

    @property
    def velocity(self):
        return self._velocity.value

    @velocity.setter
    @attribute_setter
    def velocity(self, velocity):
        self._velocity.value = velocity

    @property
    def fpa(self):
        return self._fpa.value

    @fpa.setter
    @attribute_setter
    def fpa(self, fpa):
        self._fpa.value = fpa

    @property
    def M(self):
        return self._M.value

    @M.setter
    @attribute_setter
    def M(self, M):
        self._M.value = M

    @property
    def E(self):
        return self._E.value

    @E.setter
    @attribute_setter
    def E(self, E):
        self._E.value = E

    @property
    def t_since_rp(self):
        return self._t_since_rp.value

    @t_since_rp.setter
    @attribute_setter
    def t_since_rp(self, t_since_rp):
        self._t_since_rp.value = t_since_rp

    @property
    def ascending_sign(self):
        return 1 if self.is_ascending() else -1

    @property
    def ascending(self):
        return self._is_ascending

    @ascending.setter
    @attribute_setter
    def ascending(self, ascending):
        self._is_ascending = ascending

    def is_ascending(self):

        angle_to_check = None

        if self._ta.evaluated:
            angle_to_check = self.ta
        elif self._fpa.evaluated:
            angle_to_check = self.fpa
        elif self._E.evaluated:
            angle_to_check = self.E
        elif self._velocity.evaluated:
            self._is_ascending = self.velocity.rotating().value[0].m > 0
        elif self.ascending is None:
            raise ParameterUnavailableError("Need either ta, fpa or E to check ascending")

        if angle_to_check is not None:
            self._is_ascending = 0 <= angle_to_check <= np.pi

        return self.ascending

    def angle_check_tan(self, tan_val):
        if (self.is_ascending() and tan_val >= 0) or (not self.is_ascending() and tan_val < 0):
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
        super().__init__(units.radian)
        self.orbit_requirements = [
            ('e', 'p', 'r', 'fpa'),
            ('e', 'E'),
            ('r', 'v', 'fpa'),
            ('angular_momentum', 'position', 'e_vec', 'e', 'inclination'),
            ('angular_momentum', 'position', 'ascending_node_vec', 'e', 'inclination'),
            ('position', 'e', 'inclination'),
            ('angular_momentum', 'position', 'e_vec', 'e', 'inclination', 'ascending_node_vec'),
        ]

    @orbit_setter
    def set(self, state, orbit):

        # acos((1/e) (p/r - 1))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            cos_val = (1. / orbit.e) * (orbit.p / state.r - 1.)
            self.value = state.ascending_sign * np.arccos(cos_val)

        # acos(2 atan(sqrt( (1+e)/(1-e) tan(E/2) )))
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            ta = 2 * np.arctan(np.sqrt((1 + orbit.e) / (1 - orbit.e)) * np.tan(state.E / 2))
            self.value = state.angle_check_tan(ta)

        # atan( ((rv2)/mu cos(g)sin(g)) / ((rv2)/mu cos^2(g)-1) )
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            t = (state.r * state.v ** 2) / orbit.central_body.mu
            ta = np.arctan((t * np.cos(state.fpa) * np.sin(state.fpa)) /
                           (t * np.cos(state.fpa) ** 2 - 1))
            self.value = state.angle_check_tan(ta)

        elif self.satisfied(state, orbit, self.orbit_requirements[3]) and orbit.equitorial() and not orbit.circular():

            h = orbit.angular_momentum.inertial().unit()
            e = orbit.e_vec
            r = state.position.inertial()
            self.value = np.arctan2(h.dot(e.cross(r)), r.dot(e))

        elif self.satisfied(state, orbit, self.orbit_requirements[4]) and not orbit.equitorial() and orbit.circular():
            h = orbit.angular_momentum.inertial().unit()
            r = state.position.inertial()
            n = orbit.ascending_node_vec.inertial()
            self.value = np.arctan2(r.dot(np.cross(h, n)), r.dot(n))

        elif self.satisfied(state, orbit, self.orbit_requirements[5]) and orbit.equitorial() and orbit.circular():
            r = state.position.inertial()
            self.value = np.arctan2(r[1], r[0])

        elif self.satisfied(state, orbit, self.orbit_requirements[6]) and not orbit.equitorial() and not orbit.circular():

            # unit vectors
            ehat = orbit.e_vec.inertial().unit()
            hhat = orbit.angular_momentum.inertial().unit()
            phat = np.cross(hhat, ehat)

            r = state.position.inertial().unit()
            y = r.dot(phat)
            x = r.dot(ehat)

            self.value = np.arctan2(y, x)

class PositionMagnitude(StateValue):
    symbol = 'r'

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('p', 'e', 'ta'),
            ('a', 'e', 'E'),
            ('a', 'E', 'H'),
            ('position')
        ]

    @orbit_setter
    def set(self, state, orbit):

        # p/1+ecos(ta)
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1.0 + orbit.e * np.cos(state.ta))

        # a (1-ecosE)
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 - orbit.e * np.cos(state.E))

        # |a|(E cosh H - 1)
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            self.value = abs(orbit.a) * (state.E * np.cosh(state.H) - 1)

        # |r|
        elif self.satisfied(state, orbit, self.orbit_requirements[3]):
            self.value = np.linalg.norm(state.position.value)


class ArgumentOfLatitude(StateValue):
    symbol = 'arg_latitude'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('arg_periapsis', 'ta'),
            ('inclination', 'position', 'angular_momentum')
        ]

    @orbit_setter
    def set(self, state, orbit):
        # omega + ta
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.arg_periapsis + state.ta

        # from DCM
        if self.satisfied(state, orbit, self.orbit_requirements[1]):
            rhat = state.position.inertial().unit()
            hhat = orbit.angular_momentum.inertial().unit()
            theta_hat = np.cross(hhat, rhat)

            # special case i=0
            if orbit.inclination == 0 * units.rad:
                self.value = 0
            else:
                sin_val = np.arcsin(rhat[2] / np.sin(orbit.inclination.to(units.rad)))
                cos_val = np.arccos(theta_hat[2] / np.sin(orbit.inclination.to(units.rad)))
                self.value = conics_utils.common_val(sin_val, cos_val)


class FlightPathAngle(StateValue):
    symbol = 'fpa'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('h', 'r', 'v'),
            ('velocity')
        ]

    @orbit_setter
    def set(self, state, orbit):

        # acos(h/(rv))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = state.ascending_sign * np.arccos(orbit.h / (state.r * state.v))

        # atan(vr/vt)
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            v = state.velocity.orbit_fixed()
            fpa = np.arctan(v[0] / v[1])
            self.value = state.angle_check_tan(fpa)

        if self.value and np.isnan(self.value):
            self.value = 0

class VelocityMagnitude(StateValue):
    symbol = 'v'

    def __init__(self):
        super().__init__(units.km / units.second)
        self.orbit_requirements = [
            ('r', 'se'),
            ('r', 'a'),
            ('velocity',)
        ]

    @orbit_setter
    def set(self, state, orbit):

        # sqrt( 2(se + mu/r))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(2 * (orbit.se + orbit.central_body.mu / state.r))
        # sqrt( mu(2/r-1/a) )
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = np.sqrt(orbit.central_body.mu * (2. / state.r - 1 / orbit.a))
        # |v|
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            self.value = state.velocity.norm()


class PositionVector(StateValue):
    symbol = 'position'

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('r')
        ]

    @orbit_setter
    def set(self, state, orbit):
        # r (rhat)
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = frames.Vector(orbit,
                                       state,
                                       np.array([state.r.m, 0, 0]),
                                       frames.RotatingFrame)


class VelocityVector(StateValue):
    symbol = 'velocity'

    def __init__(self):
        super().__init__(units.km / units.second)
        self.orbit_requirements = [
            ('ta', 'e', 'h'),
            ('fpa' 'v')
        ]

    @orbit_setter
    def set(self, state, orbit):

        # Rotating Frame
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            vr = (orbit.central_body.mu * orbit.e) / orbit.h * np.sin(state.ta)
            vt = orbit.central_body.mu / orbit.h * (1 + orbit.e * np.cos(state.ta))
            self.value = frames.Vector(orbit,
                                       state,
                                       np.array([vr.m, vt.m, 0]) * vr.u,
                                       frames.RotatingFrame)

        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = frames.Vector(orbit,
                                       state,
                                       np.array([state.v * np.sin(state.fpa), state.v * np.cos(state.fpa), 0]),
                                       frames.RotatingFrame)


class TimeSincePeriapsis(StateValue):
    symbol = 't_since_rp'

    def __init__(self):
        super().__init__(units.seconds)
        self.orbit_requirements = [
            ('E', 'e', 'n'),
        ]

    @orbit_setter
    def set(self, state, orbit):
        # Rotating Frame
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = (state.E - orbit.e * np.sin(state.E)) / orbit.n


class MeanAnomaly(StateValue):
    symbol = 'M'

    def __init__(self):
        super().__init__(units.radians)
        self.orbit_requirements = [
            ('n', 't_since_rp'),
            ('E', 'e')
        ]

    @orbit_setter
    def set(self, state, orbit):
        # n * ttp
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = orbit.n * state.t_since_rp

        # E - e sin(E)
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = state.E - orbit.e * np.sin(state.E)


class EccentricAnomaly(StateValue):
    symbol = 'E'

    def __init__(self):
        super().__init__(units.radians)
        self.orbit_requirements = [
            ('a', 'r', 'e'),
            ('e', 'M')
        ]

    @orbit_setter
    def set(self, state, orbit):
        # acos((a-r)/(ae))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            cos_val = (orbit.a - state.r) / (orbit.a * orbit.e)
            if cos_val > 1 and np.isclose(cos_val, 1):
                self.value = 0
            else:
                self.value = state.ascending_sign * np.arccos(
                    (orbit.a - state.r) / (orbit.a * orbit.e))

        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = self._iterative_eccentric_anomaly(state, orbit)

    def _iterative_eccentric_anomaly(self, state, orbit, **kwargs):
        return self._find_eccentric_anomaly_recursively(orbit, state.M, state.M, **kwargs)

    def _find_eccentric_anomaly_recursively(self, orbit, E, M, tol=1e-12):
        E1 = E - (E - orbit.e * np.sin(E) - M) / (1 - orbit.e * np.cos(E))
        dE = np.fabs(E - E1)

        if dE < tol:
            return E

        return self._find_eccentric_anomaly_recursively(orbit, E1, M, tol=tol)
