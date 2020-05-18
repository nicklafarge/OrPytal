from orpytal import frames, units
import numpy as np


def compare_value(orbit_or_state, poliastro_orbit, symbol, psymbol=None):
    my_value = getattr(orbit_or_state, symbol)

    if isinstance(my_value, frames.Vector):
        my_value = my_value.inertial(orbit_or_state).value

    other_value = getattr(poliastro_orbit, symbol if not psymbol else psymbol)
    unit_str = my_value.units.format_babel()

    if unit_str != 'dimensionless':
        other_value = other_value.to(unit_str)

    if my_value.units == units.rad:
        from orpytal.utils.conics_utils import angle_positive
        my_value = angle_positive(my_value)
        other_value = angle_positive(other_value.value * units.rad)

    if hasattr(other_value, 'm'):
        other_value_magnitude = other_value.m
    else:
        other_value_magnitude = other_value.value

    close_check = np.isclose(my_value.to(unit_str).m, other_value_magnitude)
    try:
        same = bool(close_check)
    except ValueError:
        same = all(close_check)

    return same
