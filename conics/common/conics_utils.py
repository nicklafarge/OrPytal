import numpy as np


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
