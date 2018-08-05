import logging
import numpy as np

from errors import ParameterUnavailableError

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

        if value_before != value_after:
            logging.debug('Set {} to {}'.format(orbit_value.symbol, orbit_value.value))

        # if value_after and isinstance(value_after, float) and np.isnan(value_after.m):
        #     logging.warning('Something went wrong setting {} (setting this value to zero now).'.format(orbit_value.symbol))
        #     orbit_value.value = 0

        # Return true if the value of the parameter has changed as as result of the function call
        return value_before != value_after

    return wrapper


def check_satisfied(obj, req):
    return hasattr(obj, '_' + req) and getattr(obj, '_' + req).evaluated

def state_orbit_satisfied(state, orbit, requirements):
    if isinstance(requirements, str):
        return check_satisfied(state, requirements) or check_satisfied(orbit, requirements)
    else:
        return all([check_satisfied(state, req) or check_satisfied(orbit, req) for req in requirements])