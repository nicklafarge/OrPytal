########### Standard ###########

########### Local ###########
from orpytal.common import units
from orpytal.errors import InvalidInputError
from orpytal.utils import conics_utils


########### External ###########


class OrbitBase(object):
    symbol = ''

    def __init__(self, units):
        self.units = units

        self.requirements = None
        self._value = None
        self.evaluated = False

    def __str__(self):
        return '{}: {}'.format(self.symbol, self._value)

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value=None):
        temp_value = conics_utils.add_units(value, self.units)

        # Resolve signs for units
        if self.units == units.rad:
            temp_value = conics_utils.angle_positive(temp_value)

        # Validate input
        if self.validate_input(temp_value):
            self._value = temp_value
            self.evaluated = True
        else:
            raise InvalidInputError()

    def __eq__(self, other):
        if isinstance(other, OrbitBase):
            return self._value == other._value
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, OrbitBase):
            return self._value != other._value
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, OrbitBase):
            return self._value > other._value
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, OrbitBase):
            return self._value >= other._value
        return NotImplemented

    def _lt__(self, other):
        if isinstance(other, OrbitBase):
            return self._value < other._value
        return NotImplemented

    def _e__(self, other):
        if isinstance(other, OrbitBase):
            return self._value <= other._value
        return NotImplemented

    def check_satisfied(self, obj, req):
        return conics_utils.check_satisfied(obj, req)

    def state_orbit_satisfied(self, state, requirements):
        return conics_utils.state_orbit_satisfied(state, requirements)

    def validate_input(self, value):
        return True
