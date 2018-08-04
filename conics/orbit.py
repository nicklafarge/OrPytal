########### Standard ###########
import logging

########### Local ###########
from base import OrbitBase
from common import units, Q_, orbit_setter
import conics_utils
import frames
from state import KeplarianState
import planet_constants
from trajectory import Trajectory

########### External ###########
import numpy as np
import scipy as sp


class Orbit():
    def __init__(self, central_body, name='', **kwargs):
        if isinstance(central_body, str):
            self.central_body = planet_constants.BODIES[central_body.upper()]
        else:
            self.central_body = central_body
        self.name = name

        self._p = SemiLatusRectum()
        self._e = Eccentricity()
        self._h = AngularMomentumMagnitude()
        self._angular_momentum = AngularMomentumVector()

        # Orientation angles
        self._arg_periapsis = ArgumentOfPeriapsis()
        self._ascending_node = LongitudeOfAscendingNode()
        self._inclination = Inclination()
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
            self._ascending_node,
            self._inclination,
            self._i,
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
        x = ['%s Orbit Info' % self.name]
        for var in self.vars:
            if var.evaluated:
                x.append(str(var))

        return '\n'.join(x)

    def propagate_full_orbit(self, state, step=0.1):
        ta_range = list(np.arange(state.ta, state.ta + 2 * np.pi, step))
        ta_range.append(2*np.pi)

        st_list = []

        for ta in ta_range:
            st = KeplarianState(self)
            st._ta.value = ta * units.radians
            st._r.set(st, self)
            st._position.set(st, self)
            st._arg_latitude.set(st, self)
            st._t_since_rp.set(st, self)
            st_list.append(st)

        return Trajectory(st_list)


    @property
    def a(self):
        return self._a.value

    @a.setter
    def a(self, a):
        self._a.value = a
        self.set_vars()

    @property
    def b(self):
        return self._b.value

    @b.setter
    def b(self, b):
        self._b.value = b
        self.set_vars()

    @property
    def period(self):
        return self._period.value

    @period.setter
    def period(self, period):
        self._period.value = period
        self.set_vars()

    @property
    def n(self):
        return self._n.value

    @n.setter
    def n(self, n):
        self._n.value = n
        self.set_vars()

    @property
    def se(self):
        return self._se.value

    @se.setter
    def se(self, se):
        self._se.value = se
        self.set_vars()

    @property
    def rp(self):
        return self._rp.value

    @rp.setter
    def rp(self, rp):
        self._rp.value = rp
        self.set_vars()

    @property
    def ra(self):
        return self._ra.value

    @ra.setter
    def ra(self, ra):
        self._ra.value = ra
        self.set_vars()

    @property
    def e(self):
        return self._e.value

    @e.setter
    def e(self, e):
        self._e.value = e

    @property
    def p(self):
        return self._p.value

    @p.setter
    def p(self, p):
        self._p.value = p

    @property
    def h(self):
        return self._h.value

    @h.setter
    def h(self, h):
        self._h.value = h

    @property
    def angular_momentum(self):
        return self._angular_momentum.value

    @angular_momentum.setter
    def angular_momentum(self, angular_momentum):
        self._angular_momentum.value = angular_momentum
        self.set_vars()

    @property
    def arg_periapsis(self):
        return self._arg_periapsis.value

    @arg_periapsis.setter
    def arg_periapsis(self, arg_periapsis):
        self._arg_periapsis.value = arg_periapsis
        self.set_vars()

    @property
    def ascending_node(self):
        return self._ascending_node.value

    @ascending_node.setter
    def ascending_node(self, ascending_node):
        self._ascending_node.value = ascending_node
        self.set_vars()

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
    def inclination(self, inclination):
        self._inclination.value = inclination
        self.set_vars()

    def set_vars(self):
        for var in self.vars:
            new_value_set = var.set(self)
            if new_value_set:
                self.set_vars()
                break

    def from_state(self, state):
        for var in self.vars:
            new_value_set = var.set_from_state(state, self)
            if new_value_set:
                self.set_vars()
                self.from_state(state)
                return True


class OrbitValue(OrbitBase):
    def set(self, orbit):
        raise NotImplementedError()

    def set_from_state(self, state, orbit):
        return False

    def satisfied(self, orbit, requirements):
        return all([self.check_satisfied(orbit, req) for req in requirements])


class SemiLatusRectum(OrbitValue):
    symbol = 'p'

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('a', 'e'),
            ('h')
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

    @orbit_setter
    def set_from_state(self, state, orbit):

        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.r * (1 + orbit.e * np.cos(state.ta))


class Eccentricity(OrbitValue):
    symbol = 'e'

    def __init__(self):
        super().__init__(units.dimensionless)
        self.orbit_requirements = [
            ('se', 'h'),
            ('a', 'rp'),
        ]
        self.orbit_state_requirements = [
            ('r', 'v', 'fpa')
        ]

    @orbit_setter
    def set(self, orbit):

        # sqrt(1 + (2 se h^2) / (mu ^ 2) )
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(1 + (2 * orbit.se * orbit.h ** 2) / (orbit.central_body.mu ** 2))

        # 1 - rp/a
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = 1.0 - orbit.rp / orbit.a

    @orbit_setter
    def set_from_state(self, state, orbit):

        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.sqrt(((state.r * state.v ** 2) / orbit.central_body.mu - 1) ** 2 * np.cos(state.fpa) ** 2 + \
                                 np.sin(state.fpa) ** 2)


class AngularMomentumMagnitude(OrbitValue):
    symbol = 'h'

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
        elif self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.linalg.norm(orbit.angular_momentum)

    @orbit_setter
    def set_from_state(self, state, orbit):

        # | r cross v |
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.linalg.norm(np.cross(state.position.value, state.velocity.value))


class AngularMomentumVector(OrbitValue):
    symbol = 'angular_momentum'

    def __init__(self):
        super().__init__(units.km ** 2 / units.s)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
            ('position', 'velocity', 'arg_periapsis', 'inclination', 'ascending_node')
        ]

    @orbit_setter
    def set(self, orbit):
        return False

    @orbit_setter
    def set_from_state(self, state, orbit):
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            h = np.cross(state.position.inertial(), state.velocity.inertial())
            self.value = frames.Vector(orbit, state, h, frames.InertialFrame)


class ArgumentOfPeriapsis(OrbitValue):
    symbol = 'arg_periapsis'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [

        ]
        self.orbit_state_requirements = [
            ('arg_latitude', 'ta')
        ]

    @orbit_setter
    def set(self, orbit):
        return False

    @orbit_setter
    def set_from_state(self, state, orbit):
        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.arg_latitude - state.ta


class LongitudeOfAscendingNode(OrbitValue):
    symbol = 'ascending_node'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
            'angular_momentum', 'inclination'
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        if self.satisfied(orbit, self.orbit_requirements[0]):
            h_xyz = orbit.angular_momentum / np.linalg.norm(orbit.angular_momentum)
            sin_val = np.arcsin(h_xyz[0] / np.sin(orbit.inclination))
            cos_val = np.arccos(-h_xyz[1] / np.sin(orbit.inclination))
            self.value = conics_utils.common_val(sin_val, cos_val)


class Inclination(OrbitValue):
    symbol = 'inclination'

    def __init__(self):
        super().__init__(units.rad)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        return False


class SemimajorAxis(OrbitValue):
    symbol = 'a'

    def __init__(self):
        super().__init__(units.km)
        self.orbit_requirements = [
            ('ra', 'rp'),
            ('se'),
            ('p', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 0.5 * (orbit.ra + orbit.rp)

        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = -orbit.central_body.mu / (2. * orbit.se)

        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = orbit.p / (1. - orbit.e ** 2)

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class SemiminorAxis(OrbitValue):
    symbol = 'b'

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

    def __init__(self):
        super().__init__(units.seconds)
        self.orbit_requirements = [
            ('n')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 2.0 * np.pi / orbit.n

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class MeanMotion(OrbitValue):
    symbol = 'n'

    def __init__(self):
        super().__init__(units.rad / units.second)
        self.orbit_requirements = [
            ('a')
        ]
        self.orbit_state_requirements = [
        ]

    @orbit_setter
    def set(self, orbit):
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.central_body.mu / (orbit.a ** 3))

    @orbit_setter
    def set_from_state(self, state, orbit):
        return False


class SpecificEnergy(OrbitValue):
    symbol = 'se'

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
