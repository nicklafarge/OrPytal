########### Standard ###########
import logging
from datetime import datetime

########### Local ###########
from orpytal import frames, output, OrbitType
from orpytal.base import OrbitBase
from orpytal.common import units
from orpytal.errors import ParameterUnavailableError, InvalidInputError
from orpytal.utils.conics_utils import *

########### External ###########
import numpy as np
from scipy.optimize import newton


class KeplarianState(object):
    def __init__(self, orbit, name='', **kwargs):
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
        self._H = HyperbolicAnomaly()

        self._ascending = None

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
            self._E,
            self._H
        ]

        if "ascending" in kwargs:
            self.ascending = kwargs.pop("ascending")
        else:
            self._ascending = None

        for k, v in kwargs.items():
            try:
                attribute = getattr(self, '_' + k)
                attribute.value = v
            except AttributeError as e:
                logging.error('KeplarianState has no parameter named {}.'.format(k))
                raise e

        self.set_vars()

    def __str__(self):
        return output.output_state(self)

    def set_vars(self):
        self.update_ascending()

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
    def H(self):
        return self._H.value

    @H.setter
    @attribute_setter
    def H(self, H):
        self._H.value = H

    @property
    def t_since_rp(self):
        return self._t_since_rp.value

    @t_since_rp.setter
    @attribute_setter
    def t_since_rp(self, t_since_rp):
        self._t_since_rp.value = t_since_rp
        if check_satisfied(self.orbit, "period"):
            self._t_since_rp.value = self._t_since_rp.value % self.orbit.period

    @property
    def ascending_sign(self):
        if self.ascending is None:
            raise ParameterUnavailableError()

        return 1 if self.ascending else -1

    @property
    def ascending(self):
        return self._ascending

    @ascending.setter
    @attribute_setter
    def ascending(self, ascending):
        if isinstance(ascending, bool):
            self._ascending = ascending
        else:
            raise InvalidInputError("Ascending must be True/False")

    def update_ascending(self):

        angle_to_check = None

        if self._ta.evaluated:
            angle_to_check = self.ta
        elif self._fpa.evaluated:
            angle_to_check = self.fpa
        elif self._E.evaluated:
            angle_to_check = self.E
        elif self._H.evaluated:
            angle_to_check = self.H
        elif self._velocity.evaluated:
            try:
                self._ascending = self.velocity.rotating().value[0].m > 0
            except ParameterUnavailableError as e:
                pass
        elif self._position.evaluated:
            try:
                self._ascending = self.position.orbit_fixed().value[1].m > 0
            except ParameterUnavailableError as e:
                pass
        elif self.ascending is None:
            logging.debug("Need either ta, fpa or E to check ascending")

        if angle_to_check is not None:
            self._ascending = 0 <= angle_to_check <= np.pi

        return self.ascending

    def angle_check_tan(self, tan_val):
        if (self.ascending and tan_val >= 0) or (not self.ascending and tan_val < 0):
            return tan_val
        else:
            return np.pi + tan_val


class StateValue(OrbitBase):
    def set(self, state, orbit):
        raise NotImplementedError()

    def satisfied(self, state, orbit, requirements):
        return self.state_orbit_satisfied(state, orbit, requirements)

    def validate_state_input(self, value, state):
        pass  # No op


class TrueAnomaly(StateValue):
    symbol = 'ta'
    name = 'True Anomaly'

    def __init__(self):
        super().__init__(units.radian)
        self.orbit_requirements = [
            ('e', 'p', 'r', 'fpa'),
            ('e', 'E'),
            ('r', 'v', 'fpa'),
            ('angular_momentum', 'position', 'e_vec', 'e', 'inclination'),
            ('angular_momentum', 'position', 'raan_vec', 'e', 'inclination'),
            ('position', 'e', 'inclination'),
            ('angular_momentum', 'position', 'e_vec', 'e', 'inclination', 'raan_vec'),
        ]

    @orbit_setter
    def set(self, state, orbit):

        if orbit.circular():
            self.value = np.nan

        # acos((1/e) (p/r - 1))
        elif self.satisfied(state, orbit, self.orbit_requirements[0]):
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
            n = orbit.raan_vec.inertial()
            self.value = np.arctan2(r.dot(np.cross(h, n)), r.dot(n))

        elif self.satisfied(state, orbit, self.orbit_requirements[5]) and orbit.equitorial() and orbit.circular():
            r = state.position.inertial()
            self.value = np.arctan2(r[1], r[0])

        elif self.satisfied(state, orbit,
                            self.orbit_requirements[6]) and not orbit.equitorial() and not orbit.circular():

            # unit vectors
            ehat = orbit.e_vec.inertial().unit()
            hhat = orbit.angular_momentum.inertial().unit()
            phat = np.cross(hhat, ehat)

            r = state.position.inertial().unit()
            y = r.dot(phat)
            x = r.dot(ehat)

            self.value = np.arctan2(y, x)

    def validate_state_input(self, value, state):
        if state.ascending is not None:
            if state.ascending and 2 * np.pi > value > np.pi:
                raise InvalidInputError("pi < value < 2pi (should be ascending orbit)")
            elif not state.ascending and np.pi > value > 0:
                raise InvalidInputError("0 < value < pi (should be descending orbit)")


class PositionMagnitude(StateValue):
    symbol = 'r'
    name = 'Pos. Mag. |r|'

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('p', 'e', 'ta'),
            ('a', 'e', 'E'),
            ('a', 'E', 'H'),
            ('position'),
            ('v', 'fpa', 'ta')
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

    def validate_state_input(self, value, state):
        if self.satisfied(state, state.orbit, ["rp"]) and value < state.orbit.rp:
            raise InvalidInputError("r < rp")

        if self.satisfied(state, state.orbit, ["ra"]) and value > state.orbit.ra:
            raise InvalidInputError("r > ra")


class ArgumentOfLatitude(StateValue):
    symbol = 'arg_latitude'
    name = 'Arg. of Latitude'

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
                self.value = common_val(sin_val, cos_val)


class FlightPathAngle(StateValue):
    symbol = 'fpa'
    name = 'Flight Path Angle'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('h', 'r', 'v'),
            ('velocity')
        ]

    @orbit_setter
    def set(self, state, orbit):
        if orbit.circular():
            self.value = 0.0 * units.rad

        # acos(h/(rv))
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = state.ascending_sign * np.arccos(orbit.h / (state.r * state.v))

        # atan(vr/vt)
        if self.satisfied(state, orbit, self.orbit_requirements[1]):
            v = state.velocity.orbit_fixed()
            fpa = np.arctan(v[0] / v[1])
            self.value = state.angle_check_tan(fpa)

        # Convert "NaN" into 0
        if self.value and np.isnan(self.value):
            self.value = 0


class VelocityMagnitude(StateValue):
    symbol = 'v'
    name = 'Vel. Mag. |v|'

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

    def validate_state_input(self, value, state):
        if self.satisfied(state, state.orbit, ["rp"]):
            v_max = state.orbit.get_state(r=state.orbit.rp).v
            if value > v_max:
                raise InvalidInputError("v > v_max (%0.6e)".format(v_max))

        if self.satisfied(state, state.orbit, ["ra"]):
            v_min = state.orbit.get_state(r=state.orbit.ra).v
            if value < v_min:
                raise InvalidInputError("v < v_min (%0.6e)".format(v_min))


class PositionVector(StateValue):
    symbol = 'position'
    name = 'Position'

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
    name = 'Velocity'

    def __init__(self):
        super().__init__(units.km / units.second)
        self.orbit_requirements = [
            ('ta', 'e', 'h'),
            ('fpa', 'v')
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
                                       np.array([state.v.m * np.sin(state.fpa), state.v.m * np.cos(state.fpa), 0]),
                                       frames.RotatingFrame)


class TimeSincePeriapsis(StateValue):
    symbol = 't_since_rp'
    name = 'Time Since Peri.'

    def __init__(self):
        super().__init__(units.seconds)
        self.orbit_requirements = [
            ('E', 'e', 'n'),
        ]

    @orbit_setter
    def set(self, state, orbit):
        # ( E - e sin(E) ) / n
        if self.satisfied(state, orbit, self.orbit_requirements[0]):
            self.value = (state.E - orbit.e * np.sin(state.E)) / orbit.n


class MeanAnomaly(StateValue):
    symbol = 'M'
    name = 'Mean Anomaly'

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
    name = 'Eccentric Anomaly'

    def __init__(self):
        super().__init__(units.radians)
        self.orbit_requirements = [
            ('a', 'r', 'e'),
            ('e', 'M'),
            ('ta', 'e')
        ]

    @orbit_setter
    def set(self, state, orbit):
        if orbit.circular():
            self.value = np.nan

        # Only proceed if elliptic orbit
        if orbit.type() != OrbitType.Elliptic:
            pass

        # acos((a-r)/(ae))
        elif self.satisfied(state, orbit, self.orbit_requirements[0]):
            cos_val = (orbit.a - state.r) / (orbit.a * orbit.e)
            if cos_val > 1 and np.isclose(cos_val, 1):
                self.value = 0
            else:
                self.value = state.ascending_sign * np.arccos(
                    (orbit.a - state.r) / (orbit.a * orbit.e))

        # Newton raphson to find eccentric anomaly from from mean anomaly
        elif self.satisfied(state, orbit, self.orbit_requirements[1]):
            self.value = newton(lambda E, state: E - state.orbit.e.m * np.sin(E) - state.M.m,
                                state.M.m,
                                args=(state,),
                                fprime=lambda E, state: 1 - state.orbit.e.m * np.cos(E),
                                tol=1e-12)

        # acos(2 atan(sqrt( (1+e)/(1-e) tan(E/2) )))
        elif self.satisfied(state, orbit, self.orbit_requirements[2]):
            E = 2 * np.arctan(np.tan(state.ta / 2) / np.sqrt((1 + orbit.e) / (1 - orbit.e)))
            self.value = state.angle_check_tan(E)


class HyperbolicAnomaly(StateValue):
    symbol = 'H'
    name = 'Hyperbolic Anomaly'

    def __init__(self):
        super().__init__(units.radians)
        self.orbit_requirements = [
            ('e', 'a', 'n')
        ]

    @orbit_setter
    def set(self, state, orbit):
        if orbit.type() != OrbitType.Hyperbolic:
            pass

        # TODO: fix this!
        # Newton raphson to find hyperbolic anomaly
        # elif self.satisfied(state, orbit, self.orbit_requirements[0]):
        #     def findH(H, state):
        #         e = state.orbit.e
        #         mu = state.orbit.central_body.mu
        #         a = state.orbit.a
        #         n = state.orbit.n
        #         ttp = 1 / n * (e * np.sinh(H) - H)
        #         return e * np.sinh(H) - H - np.sqrt(mu / (np.abs(a) ** 3)) * ttp
        #
        #     self.value = newton(findH,
        #                         0.1,  # TODO is there a better initial guess?
        #                         args=(state,),
        #                         tol=1e-12)
