########### Standard ###########
import logging

########### Local ###########
from orpytal.base import OrbitBase
from orpytal.common import units, orbit_setter, attribute_setter
from orpytal.errors import ParameterUnavailableError
from orpytal import frames, output
from orpytal.state import KeplarianState
from orpytal.trajectory import Trajectory

########### External ###########
import numpy as np


class Orbit(object):
    def __init__(self, central_body, name='', **kwargs):
        self.central_body = central_body
        self.name = name

        self._p = SemiLatusRectum()
        self._e = Eccentricity()
        self._e_vec = EccentricityVector()
        self._h = AngularMomentumMagnitude()
        self._angular_momentum = AngularMomentumVector()

        # Orientation angles
        self._arg_periapsis = ArgumentOfPeriapsis()
        self._raan = LongitudeOfAscendingNode()
        self._inclination = Inclination()
        self._raan_vec = AscendingNodeVector()
        self._i = self._inclination

        self._a = SemimajorAxis()
        self._b = SemiminorAxis()
        self._period = OrbitalPeriod()
        self._n = MeanMotion()
        self._se = SpecificEnergy()
        self._rp = Periapsis()
        self._ra = Apoapsis()

        self.vars = [
            self._p,
            self._e,
            self._e_vec,
            self._h,
            self._angular_momentum,
            self._a,
            self._b,
            self._period,
            self._n,
            self._se,
            self._rp,
            self._ra,
            self._arg_periapsis,
            self._raan,
            self._raan_vec,
            self._inclination
        ]

        for k, v in kwargs.items():
            try:
                attribute = getattr(self, '_' + k)
                attribute.value = v
            except AttributeError as e:
                logging.error('Orbit has no parameter named {}.'.format(k))
                raise e

        self.set_vars()

    def __str__(self):
        return output.output_orbit(self)

    def propagate_full_orbit(self, n=150):
        t_range = np.linspace(start=0, stop=self.period.m, num=n)

        st_list = []
        for t in t_range:
            st = KeplarianState(self)
            st.t_since_rp = t
            st._t_since_rp.value = t
            st._M.set(st, self)
            st._E.set(st, self)
            st._r.set(st, self)
            st._ta.set(st, self)
            st._arg_latitude.set(st, self)
            st._position.set(st, self)
            st._velocity.set(st, self)

            st_list.append(st)

        return Trajectory(st_list)

    def propagate_full_orbit_ta(self, state, step=0.1):
        ta_range = list(np.arange(state.ta, state.ta + 2 * np.pi, step))
        ta_range.append(state.ta.m + 2 * np.pi)

        st_list = []

        for ta in ta_range:
            st = KeplarianState(self)
            st._ta.value = ta
            st._r.set(st, self)
            st._position.set(st, self)
            st._arg_latitude.set(st, self)
            st._t_since_rp.set(st, self)
            st_list.append(st)

        return Trajectory(st_list)

    def set_vars(self):
        for var in self.vars:
            new_value_set = var.set(self)
            if new_value_set:
                logging.debug('Set {} to {}'.format(var.symbol, var.value))
                self.set_vars()
                break

    def from_state(self, state):
        for var in self.vars:
            new_value_set = var.set_from_state(state, self)
            if new_value_set:
                self.set_vars()
                self.from_state(state)
                return True

    def equitorial(self, tol=1e-8):
        if self.inclination is not None:
            return np.abs(self.inclination) < tol
        else:
            raise ParameterUnavailableError('Need inclination to see if equitorial')

    def circular(self, tol=1e-8):
        if self.e is not None:
            return self.e < tol
        else:
            raise ParameterUnavailableError('Need eccentricity to see if circular')

    def compare(self, orbit):
        all_same = True

        # Compare which vars are evaluated
        evaluated_vars = sorted([v.symbol for v in self.vars if v.evaluated])
        other_evaluated_vars = sorted([v.symbol for v in orbit.vars if v.evaluated])
        for i, var in enumerate(evaluated_vars):
            if var != other_evaluated_vars[i]:
                logging.warning('Different Evaluated Variables!')
                return False
            # assert var == other_evaluated_vars[i]

        # Compare values of evaluated variables
        for var in self.vars:
            if var.evaluated and hasattr(orbit, var.symbol):
                if isinstance(var.value, units.Quantity):
                    same = np.isclose(var.value, getattr(orbit, var.symbol))
                elif isinstance(var.value, frames.Vector):

                    vec1 = var.value
                    vec2 = getattr(orbit, var.symbol)
                    same = vec1 == vec2
                else:
                    same = False

            if same:
                logging.debug('Checked {} [✓]'.format(var.symbol))
            else:
                all_same = False
                logging.warning('Error Found for {} [x]'.format(var.symbol))
                logging.warning('My value: {}'.format(var.value))
                logging.warning('Their value: {}'.format(getattr(orbit, var.symbol)))

        return all_same

    @property
    def a(self):
        return self._a.value

    @a.setter
    @attribute_setter
    def a(self, a):
        self._a.value = a

    @property
    def b(self):
        return self._b.value

    @b.setter
    @attribute_setter
    def b(self, b):
        self._b.value = b

    @property
    def period(self):
        return self._period.value

    @period.setter
    @attribute_setter
    def period(self, period):
        self._period.value = period

    @property
    def n(self):
        return self._n.value

    @n.setter
    @attribute_setter
    def n(self, n):
        self._n.value = n

    @property
    def se(self):
        return self._se.value

    @se.setter
    @attribute_setter
    def se(self, se):
        self._se.value = se

    @property
    def rp(self):
        return self._rp.value

    @rp.setter
    @attribute_setter
    def rp(self, rp):
        self._rp.value = rp

    @property
    def ra(self):
        return self._ra.value

    @ra.setter
    @attribute_setter
    def ra(self, ra):
        self._ra.value = ra

    @property
    def e(self):
        return self._e.value

    @e.setter
    @attribute_setter
    def e(self, e):
        self._e.value = e

    @property
    def e_vec(self):
        return self._e_vec.value

    @e_vec.setter
    @attribute_setter
    def e_vec(self, e_vec):
        self._e_vec.value = e_vec

    @property
    def p(self):
        return self._p.value

    @p.setter
    @attribute_setter
    def p(self, p):
        self._p.value = p

    @property
    def h(self):
        return self._h.value

    @h.setter
    @attribute_setter
    def h(self, h):
        self._h.value = h

    @property
    def angular_momentum(self):
        return self._angular_momentum.value

    @angular_momentum.setter
    @attribute_setter
    def angular_momentum(self, angular_momentum):
        self._angular_momentum.value = angular_momentum

    @property
    def arg_periapsis(self):
        return self._arg_periapsis.value

    @arg_periapsis.setter
    @attribute_setter
    def arg_periapsis(self, arg_periapsis):
        self._arg_periapsis.value = arg_periapsis

    @property
    def raan(self):
        return self._raan.value

    @raan.setter
    @attribute_setter
    def raan(self, raan):
        self._raan.value = raan

    @property
    def raan_vec(self):
        return self._raan_vec.value

    @raan_vec.setter
    @attribute_setter
    def raan_vec(self, raan_vec):
        self._raan_vec.value = raan_vec

    @property
    def i(self):
        return self.inclination

    @i.setter
    def i(self, i):
        self.inclination = i

    @property
    def inclination(self):
        return self._inclination.value

    @inclination.setter
    @attribute_setter
    def inclination(self, inclination):
        self._inclination.value = inclination


class OrbitValue(OrbitBase):
    def set(self, orbit):
        raise NotImplementedError()

    def set_from_state(self, state, orbit):
        return False

    def satisfied(self, orbit, requirements):
        if isinstance(requirements, str):
            return self.check_satisfied(orbit, requirements)
        else:
            return all([self.check_satisfied(orbit, req) for req in requirements])


class SemiLatusRectum(OrbitValue):
    symbol = 'p'
    name = "Semilatus Rectum"

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('a', 'e'),
            ('h'),
            ('ra', 'e')
        ]
        self.orbit_state_requirements = [
            ('r', 'e', 'ta')
        ]

    @orbit_setter
    def set(self, orbit):

        # a(1-e^2)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.a * (1 - orbit.e ** 2)

        # h^2/mu
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.h ** 2 / orbit.central_body.mu

        # ra * (1 - e)
        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = orbit.ra * (1 - orbit.e)


    @orbit_setter
    def set_from_state(self, state, orbit):

        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.r * (1 + orbit.e * np.cos(state.ta))


class Eccentricity(OrbitValue):
    symbol = 'e'
    name = "Eccentricity"

    def __init__(self):
        super().__init__(units.dimensionless)
        self.orbit_requirements = [
            ('se', 'h'),
            ('a', 'rp'),
            ('e_vec'),
            ('p', 'ra'),
            ('a', 'ra'),
            ('p', 'rp')
        ]
        self.orbit_state_requirements = [
            ('r', 'v', 'fpa')
        ]

    @orbit_setter
    def set(self, orbit):

        # sqrt(1 + (2 se h^2) / (mu ^ 2) )
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(1. + (2. * orbit.se * orbit.h ** 2.) / (orbit.central_body.mu ** 2.))

        # 1 - rp/a
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = 1. - orbit.rp / orbit.a

        # |e|
        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = np.linalg.norm(orbit.e_vec)

        # 1 - p/ra
        elif self.satisfied(orbit, self.orbit_requirements[3]):
            self.value = 1. - orbit.p / orbit.ra

        # ra/a - 1
        elif self.satisfied(orbit, self.orbit_requirements[4]):
            self.value = orbit.ra/orbit.a - 1.

        # p / rp - 1
        elif self.satisfied(orbit, self.orbit_requirements[5]):
            self.value = orbit.p / orbit.rp - 1.

    @orbit_setter
    def set_from_state(self, state, orbit):

        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.sqrt(((state.r * state.v ** 2) / orbit.central_body.mu - 1) ** 2 * np.cos(state.fpa) ** 2 + \
                                 np.sin(state.fpa) ** 2)


class EccentricityVector(OrbitValue):
    symbol = 'e_vec'
    name = "Eccentricity Vector"

    def __init__(self):
        super().__init__(units.dimensionless)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
            ('position', 'velocity')
        ]

    @orbit_setter
    def set(self, orbit):
        return

    @orbit_setter
    def set_from_state(self, state, orbit):
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            v = state.velocity.inertial()  # inertial velocity (xyz)
            r = state.position.inertial()  # inertial position (xyz)
            mu = orbit.central_body.mu  # gravitational parameter
            val = ((v.dot(v) - mu / r.norm()) * r.value - r.dot(v) * v.value) / mu
            self.value = frames.Vector(orbit, state, val, frames.InertialFrame)


class AngularMomentumMagnitude(OrbitValue):
    symbol = 'h'
    name = "Angular Momentum"

    def __init__(self):
        super().__init__(units.km ** 2 / units.s)
        self.orbit_requirements = [
            ('p'),
            ('angular_momentum')
        ]
        self.orbit_state_requirements = [
            ('position', 'velocity')
        ]

    @orbit_setter
    def set(self, orbit):

        # sqrt(p mu)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.p * orbit.central_body.mu)
        # | h |
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = np.linalg.norm(orbit.angular_momentum)

    @orbit_setter
    def set_from_state(self, state, orbit):

        # | r cross v |
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.linalg.norm(np.cross(state.position.value, state.velocity.value))


class AngularMomentumVector(OrbitValue):
    symbol = 'angular_momentum'
    name = "Angular Momentum Vector"

    def __init__(self):
        super().__init__(units.km ** 2 / units.s)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
            ('raan', 'inclination', 'h'),
            ('position', 'velocity')
        ]

    @orbit_setter
    def set(self, orbit):
        return False

    @orbit_setter
    def set_from_state(self, state, orbit):
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            h_vector = orbit.h * np.array([np.sin(orbit.raan) * np.sin(orbit.inclination),
                                           -np.cos(orbit.raan) * np.sin(orbit.inclination),
                                           np.cos(orbit.inclination)
                                           ])
            self.value = frames.Vector(orbit, state, h_vector, frames.InertialFrame)
        elif self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[1]):
            h = state.position.inertial().cross(state.velocity.inertial())
            self.value = frames.Vector(orbit, state, h, frames.InertialFrame)


class ArgumentOfPeriapsis(OrbitValue):
    symbol = 'arg_periapsis'
    name = "Arg of Periapsis"

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('e'),
            ('e_vec', 'e'),
            ('e_vec', 'angular_momentum', 'raan_vec')

        ]
        self.orbit_state_requirements = [
            ('arg_latitude', 'ta'),
        ]

    @orbit_setter
    def set(self, orbit):
        if self.satisfied(orbit, self.orbit_requirements[0]) and orbit.circular():
            self.value = 0

        elif self.satisfied(orbit, self.orbit_requirements[1]) and orbit.equitorial():
            self.value = np.arctan2(orbit.e_vec[1], orbit.e_vec[0])

        elif self.satisfied(orbit, self.orbit_requirements[2]) and not orbit.equitorial() and not orbit.circular():
            e = orbit.e_vec.inertial().value.m
            h_xyz = orbit.angular_momentum.inertial().unit()
            line_of_nodes = orbit.raan_vec
            self.value = np.arctan2(e.dot(np.cross(h_xyz, line_of_nodes)), e.dot(line_of_nodes))

    @orbit_setter
    def set_from_state(self, state, orbit):
        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.arg_latitude - state.ta


class LongitudeOfAscendingNode(OrbitValue):
    symbol = 'raan'
    name = "Ascending Node"

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('inclination'),
            ('raan_vec')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        # special case if i=0
        if self.satisfied(orbit, self.orbit_requirements[0]) and orbit.equitorial():
            self.value = 0

        # Angle from line of nodes
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = np.arctan2(orbit.raan_vec[1], orbit.raan_vec[0])


class AscendingNodeVector(OrbitValue):
    symbol = 'raan_vec'
    name = "Ascending Node Vector"

    def __init__(self):
        super().__init__(units.dimensionless)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
            ('angular_momentum')
        ]

    @orbit_setter
    def set(self, orbit):
        return False

    @orbit_setter
    def set_from_state(self, state, orbit):
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            h_xyz = orbit.angular_momentum.inertial().unit()
            value = np.cross([0, 0, 1], h_xyz)
            self.value = frames.Vector(orbit, state, value, frames.InertialFrame)


class Inclination(OrbitValue):
    symbol = 'inclination'
    name = "Inclination"

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            ('angular_momentum')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        if self.satisfied(orbit, self.orbit_requirements[0]):
            h_xyz = orbit.angular_momentum.unit()
            self.value = np.arccos(h_xyz[2])


class SemimajorAxis(OrbitValue):
    symbol = 'a'
    name = "Semimajor Axis"

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('ra', 'rp'),
            ('se'),
            ('p', 'e'),
            ('rp', 'e'),
            ('n')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        # (ra+rp)/2
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 0.5 * (orbit.ra + orbit.rp)

        # -mu/(2*se)
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = -orbit.central_body.mu / (2. * orbit.se)

        # p / (1 - e^2)
        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = orbit.p / (1. - orbit.e ** 2)

        # rp / (1 - e)
        elif self.satisfied(orbit, self.orbit_requirements[3]):
            self.value = orbit.rp / (1 - orbit.e)

        # (mu / n^2) ^ (1/3)
        elif self.satisfied(orbit, self.orbit_requirements[4]):
            self.value = (orbit.central_body.mu / orbit.n ** 2) ** (1./3.)

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class SemiminorAxis(OrbitValue):
    symbol = 'b'
    name = "Semiminor Axis"

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1. + orbit.e ** 2.)

        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * np.sqrt(1 - orbit.e ** 2)

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class OrbitalPeriod(OrbitValue):
    symbol = 'period'
    name = "Orbital Period"

    def __init__(self):
        super().__init__(units.seconds)
        self.orbit_requirements = [
            ('n')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        # 2pi / n
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 2.0 * np.pi / orbit.n

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class MeanMotion(OrbitValue):
    symbol = 'n'
    name = "Mean Motion"

    def __init__(self):
        super().__init__(units.rad / units.second)
        self.orbit_requirements = [
            ('a'),
            ('period')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        # sqrt(mu / a^3)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.central_body.mu / (orbit.a ** 3))

        # 2pi / period
        if self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = 2 * np.pi / orbit.period


    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class SpecificEnergy(OrbitValue):
    symbol = 'se'
    name = "Specific Energy"

    def __init__(self):
        super().__init__(units.km ** 2 / units.s ** 2)
        self.orbit_requirements = [
            ('a')
        ]
        self.orbit_state_requirements = [
            ('r', 'v')
        ]

    @orbit_setter
    def set(self, orbit):

        # -mu/(2a)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = -orbit.central_body.mu / (2 * orbit.a)

    @orbit_setter
    def set_from_state(self, state, orbit):

        # v^2/2 - mu/r
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.v ** 2. / 2. - orbit.central_body.mu / state.r


class Apoapsis(OrbitValue):
    symbol = 'ra'
    name = "Rad. Apoapsis"

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        # p / (1 - e)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1 - orbit.e)

        # a * (1 + e)
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 + orbit.e)

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class Periapsis(OrbitValue):
    symbol = 'rp'
    name = "Rad. Periapsis"

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        # p / (1 + e)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1 + orbit.e)

        # a * (1 - e)
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 - orbit.e)

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False
