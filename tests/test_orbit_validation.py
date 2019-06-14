########### Standard ###########
import logging
import itertools
import unittest

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, integration, frames, Trajectory
from orpytal.errors import InvalidInputError

########### External ###########
import numpy as np

logging.disable(logging.CRITICAL)


class TestOrbitValidationMethods(unittest.TestCase):
    def test_validation_from_constructor(self):
        """
            Test construtor kwargs with semimajor axis
        """

        assert (Orbit(bodies.earth, a=0).a == None)
        assert (Orbit(bodies.earth, a=18000).a != None)

    def test_with_zero_semimajor_axis(self):
        """
            Test with invalid semimajor axis
        """

        orbit = Orbit(bodies.earth)
        orbit.a = 0
        assert orbit.a == None

        orbit.a = 18000
        assert orbit.a != None
