import numpy as np


def loc(v1, v2, theta):
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(theta))


def loc_angle(v1, v2, v3):
    return np.arccos((v3 ** 2 - v1 ** 2 - v2 ** 2) / -(2 * v1 * v2))

