########### Standard ###########

########### Local ###########
from base import OrbitBase
from common import ureg, Q_
import conics_utils
import frames
from planet_constants import BODIES

########### External ###########
import numpy as np
import scipy as sp


class Orbit():
    def __init__(self, central_body, name=''):
        self.central_body = BODIES[central_body.upper()]
        self.name = name

        self._p = SemiLatusRectum()
        self._e = Eccentricity()
        self._h = AngularMomentumMagnitude()
        self._angular_momentum = AngularMomentumVector()

        # Orientation angles
        self._arg_periapsis = ArgumentOfPeriapsis()
        self._ascending_node = LongitudeOfAscendingNode()
        self._inclination = Inclination()

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
        ]

    def __str__(self):
        x = ['%s Orbit Info' % self.name]
        for var in self.vars:
            if var.evaluated:
                x.append(str(var))

        return '\n'.join(x)

    @property
    def a(self):
        return self._a.value

    @a.setter
    def a(self, a=None):
        self._a.value = a
        self.set_vars()

    @property
    def b(self):
        return self._b.value

    @b.setter
    def b(self, b=None):
        self._b.value = b
        self.set_vars()

    @property
    def period(self):
        return self._period.value

    @period.setter
    def period(self, period=None):
        self._period.value = period
        self.set_vars()

    @property
    def n(self):
        return self._n.value

    @n.setter
    def n(self, n=None):
        self._n.value = n
        self.set_vars()

    @property
    def se(self):
        return self._se.value

    @se.setter
    def se(self, se=None):
        self._se.value = se
        self.set_vars()

    @property
    def rp(self):
        return self._rp.value

    @rp.setter
    def rp(self, rp=None):
        self._rp.value = rp
        self.set_vars()

    @property
    def ra(self):
        return self._ra.value

    @ra.setter
    def ra(self, ra=None):
        self._ra.value = ra
        self.set_vars()

    @property
    def e(self):
        return self._e.value

    @e.setter
    def e(self, e=None):
        self._e.value = e

    @property
    def p(self):
        return self._p.value

    @p.setter
    def p(self, p=None):
        self._p.value = p

    @property
    def h(self):
        return self._h.value

    @h.setter
    def h(self, h=None):
        self._h.value = h

    @property
    def angular_momentum(self):
        return self._angular_momentum.value

    @angular_momentum.setter
    def angular_momentum(self, angular_momentum=None):
        self._angular_momentum.value = angular_momentum
        self.set_vars()

    @property
    def arg_periapsis(self):
        return self._arg_periapsis.value

    @arg_periapsis.setter
    def arg_periapsis(self, arg_periapsis=None):
        self._arg_periapsis.value = arg_periapsis
        self.set_vars()

    @property
    def ascending_node(self):
        return self._ascending_node.value

    @ascending_node.setter
    def ascending_node(self, ascending_node=None):
        self._ascending_node.value = ascending_node
        self.set_vars()

    @property
    def inclination(self):
        return self._inclination.value

    @inclination.setter
    def inclination(self, inclination=None):
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
                break


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
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('a', 'e'),
            ('h')
        ]
        self.orbit_state_requirements = [
            ('r', 'e', 'ta')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # a(1-e^2)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.a * (1 - orbit.e ** 2)

        # h^2/mu
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.h ** 2 / orbit.central_body.mu

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.r * (1 + orbit.e * np.cos(state.ta))

        # Requirements not met
        else:
            return False

        return True


class Eccentricity(OrbitValue):
    symbol = 'e'

    def __init__(self):
        super().__init__(ureg.dimensionless)
        self.orbit_requirements = [
            ('se', 'h'),
            ('a', 'rp'),
        ]
        self.orbit_state_requirements = [
            ('r', 'v', 'fpa')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # sqrt(1 + (2 se h^2) / (mu ^ 2) )
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(1 + (2 * orbit.se * orbit.h ** 2) / (orbit.central_body.mu ** 2))

        # 1 - rp/a
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = 1.0 - orbit.rp / orbit.a

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.sqrt(((state.r * state.v ** 2) / orbit.central_body.mu - 1) ** 2 * np.cos(state.fpa) ** 2 + \
                                 np.sin(state.fpa) ** 2)

        # Requirements not met
        else:
            return False

        return True


class AngularMomentumMagnitude(OrbitValue):
    symbol = 'h'

    def __init__(self):
        super().__init__(ureg.km ** 2 / ureg.s)
        self.orbit_requirements = [
            ('p')
        ]
        self.orbit_state_requirements = [
            ('pos', 'vel')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # sqrt(p mu)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.p * orbit.central_body.mu)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        # | r cross v |
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = np.linalg.norm(np.cross(state.pos.value, state.vel.value))

        # Requirements not met
        else:
            return False

        return True


class AngularMomentumVector(OrbitValue):
    symbol = 'angular_momentum'

    def __init__(self):
        super().__init__(ureg.km ** 2 / ureg.s)
        self.orbit_requirements = [
            ('p')
        ]
        self.orbit_state_requirements = [
            ('pos', 'vel', 'arg_periapsis', 'inclination', 'ascending_node')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # sqrt(p mu)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.p * orbit.central_body.mu)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            h = np.cross(state.pos.inertial(), state.vel.inertial())
            self.value = frames.Vector(h, frames.InertialFrame)

        # Requirements not met
        else:
            return False

        return True


class ArgumentOfPeriapsis(OrbitValue):
    symbol = 'arg_periapsis'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [

        ]
        self.orbit_state_requirements = [
            ('arg_latitude', 'ta')
        ]

    def set(self, orbit):
        return False

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        # r(1+ecos(ta))
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.arg_latitude - state.ta

        # Requirements not met
        else:
            return False

        return True


class LongitudeOfAscendingNode(OrbitValue):
    symbol = 'ascending_node'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [
            'angular_momentum', 'inclination'
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        if self.satisfied(orbit, self.orbit_requirements[0]):
            h_xyz = orbit.angular_momentum / np.linalg.norm(orbit.angular_momentum)
            sin_val = np.arcsin(h_xyz[0] / np.sin(orbit.inclination))
            cos_val = np.arccos(-h_xyz[1] / np.sin(orbit.inclination))
            self.value = conics_utils.common_val(sin_val, cos_val)

        # Requirements not met
        else:
            return False

        return True


class Inclination(OrbitValue):
    symbol = 'inclination'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        return False


class SemimajorAxis(OrbitValue):
    symbol = 'a'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('ra', 'rp'),
            ('se'),
            ('p', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 0.5 * (orbit.ra + orbit.rp)

        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = -orbit.central_body.mu / (2. * orbit.se)

        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = orbit.p / (1. - orbit.e ** 2)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False


class SemiminorAxis(OrbitValue):
    symbol = 'b'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1. + orbit.e ** 2.)

        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * np.sqrt(1 - orbit.e ** 2)

        elif self.satisfied(orbit, self.orbit_requirements[2]):
            self.value = orbit.p / (1. - orbit.e ** 2)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False


class OrbitalPeriod(OrbitValue):
    symbol = 'period'

    def __init__(self):
        super().__init__(ureg.seconds)
        self.orbit_requirements = [
            ('n')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = 2.0 * np.pi / orbit.n

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False


class MeanMotion(OrbitValue):
    symbol = 'n'

    def __init__(self):
        super().__init__(ureg.rad / ureg.second)
        self.orbit_requirements = [
            ('a')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = np.sqrt(orbit.central_body.mu / (orbit.a ** 3))

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False


class SpecificEnergy(OrbitValue):
    symbol = 'se'

    def __init__(self):
        super().__init__(ureg.km ** 2 / ureg.s ** 2)
        self.orbit_requirements = [
            ('a')
        ]
        self.orbit_state_requirements = [
            ('r', 'v')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # -mu/(2a)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = -orbit.central_body.mu / (2 * orbit.a)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        if self.evaluated:
            return False

        # v^2/2 - mu/r
        if self.state_orbit_satisfied(state, orbit, self.orbit_state_requirements[0]):
            self.value = state.v ** 2. / 2. - orbit.central_body.mu / state.r

        # Requirements not met
        else:
            return False

        return True


class Apoapsis(OrbitValue):
    symbol = 'ra'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # p / (1 + e)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1 + orbit.e)

        # a * (1 + e)
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 + orbit.e)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False


class Periapsis(OrbitValue):
    symbol = 'rp'

    def __init__(self):
        super().__init__(ureg.km)
        self.orbit_requirements = [
            ('p', 'e'),
            ('a', 'e')
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        # p / (1 - e)
        if self.satisfied(orbit, self.orbit_requirements[0]):
            self.value = orbit.p / (1 - orbit.e)

        # a * (1 - e)
        elif self.satisfied(orbit, self.orbit_requirements[1]):
            self.value = orbit.a * (1 - orbit.e)

        # Requirements not met
        else:
            return False

        return True

    def set_from_state(self, state, orbit):
        return False
