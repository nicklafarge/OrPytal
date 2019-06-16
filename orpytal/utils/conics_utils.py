import logging
import numpy as np

from orpytal import frames
from orpytal.errors import ParameterUnavailableError, InvalidInputError
from orpytal.common import units


def loc(v1, v2, theta):
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(theta))


def loc_angle(v1, v2, v3):
    return np.arccos((v3 ** 2 - v1 ** 2 - v2 ** 2) / -(2 * v1 * v2))


def common_val(sin_val, cos_val, tol=1e-3):
    one = sin_val
    if one > 0:
        two = np.pi - sin_val
    else:
        two = -np.pi - sin_val

    three = cos_val
    four = -cos_val

    if np.abs(one - three) < tol or np.abs(one - four) < tol:
        val = one
    elif np.abs(two - three) < tol or np.abs(two - four) < tol:
        val = two
    else:
        raise ValueError('No common value found for %f and %f\n' % (sin_val, cos_val))

    return val


def orbit_setter(setter_function):
    def wrapper(*args):
        orbit_value = args[0]

        # If this orbital parameter already has a value, return immediately
        if orbit_value.evaluated:
            return False

        value_before = orbit_value.value
        try:
            setter_function(*args)
        except ParameterUnavailableError as e:
            logging.debug('Assertion Error: {}'.format(e))
            pass

        value_after = orbit_value.value

        try:
            value_changed = bool(value_before != value_after)
        except ValueError as ve:
            value_changed = any(value_before != value_after)

        if value_changed:
            logging.debug('Set {} to {}'.format(orbit_value.symbol, orbit_value.value))

        # Return true if the value of the parameter has changed as as result of the function call
        return value_changed

    return wrapper


def set_attribute(orbit_or_state, val, setter_function, var_name):
    cls_name = orbit_or_state.__class__.__name__

    is_keplarian_state = cls_name == 'KeplarianState'

    if is_keplarian_state:
        # Allow for tuple vector inputs
        if isinstance(val, tuple) and \
                hasattr(val[0], '__len__') and \
                (isinstance(val[1], frames.CoordinateFrame) or val[1].__bases__[0] == frames.CoordinateFrame):
            val = frames.Vector(val[0], val[1])

    # Validate input
    try:
        var = getattr(orbit_or_state, "_" + var_name)
        if var.evaluated:
            logging.error("Variable {} already has a value set - it can't be set again. Create a new instance "
                          "instead of changing an existing one".format(var_name))
            return False
        elif hasattr(var, "validate_state_input"):
            var.validate_state_input(add_units(val, var.units), orbit_or_state)
    except InvalidInputError as iie:
        logging.error("Invalid input for {}. Validation failed due to:\n {}".format(var_name, str(iie)))
        return False

    setter_function(orbit_or_state, val)
    logging.debug('Set {} to {}'.format(var_name, val))

    orbit_or_state.set_vars()
    return True


def attribute_setter(setter_function):
    def wrapper(*args):
        set_attribute(*args, setter_function, setter_function.__name__)

    return wrapper


def check_satisfied(obj, req):
    return hasattr(obj, '_' + req) and getattr(obj, '_' + req).evaluated


def state_orbit_satisfied(state, requirements):
    if isinstance(requirements, str):
        return check_satisfied(state, requirements) or check_satisfied(state.orbit, requirements)
    else:
        return all([check_satisfied(state, req) or check_satisfied(state.orbit, req) for req in requirements])


def add_units(value, value_units):
    if isinstance(value, units.Quantity):
        return value.to(value_units)
    elif isinstance(value, frames.Vector):
        if isinstance(value.value, units.Quantity):
            value.value = value.value.to(value_units)
            return value
        else:
            value.value = value.value * value_units
            return value
    else:
        return value * value_units


def angle_positive(value):
    return ((value.to(units.rad) + (2 * np.pi)) % (2 * np.pi)).to('rad')


def angle_pos_neg(value):
    angle = angle_positive(value)
    if angle > np.pi:
        return angle - 2 * np.pi
    else:
        return angle
