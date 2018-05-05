########### Standard ###########

########### Local ###########
from base import OrbitBase
from common import ureg, Q_
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
        self._h = AngularMomentum()

        # Orientation angles
        self._arg_periapsis = ArgumentOfPeriapsis()
        self._ascending_node = LongitudeOfAscendingNode()
        self._inclination = Inclination()

        self.vars = [
            self._p,
            self._e,
            self._h,
            self._arg_periapsis,
            self._ascending_node,
            self._inclination,
        ]

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
            ('se', 'h')
        ]

    def set(self, orbit):
        if self.evaluated:
            return False


class AngularMomentum(OrbitValue):
    symbol = 'h'

    def __init__(self):
        super().__init__(ureg.km ** 2 / ureg.s)
        self.orbit_requirements = [
            ('p')
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
        if self.evaluated:
            return False

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
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        return False


class Inclination(OrbitValue):
    symbol = 'inclination'

    def __init__(self):
        super().__init__(ureg.rad)
        self.orbit_requirements = [
        ]
        self.orbit_state_requirements = [
        ]

    def set(self, orbit):
        if self.evaluated:
            return False

        return False
