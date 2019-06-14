########### Standard ###########
import logging
import unittest

########### Local ###########
from orpytal import Orbit, bodies, OrbitType
from orpytal.common import units


########### External ###########
import numpy as np


class TestOrbitTypeDetection(unittest.TestCase):
    def test_from_eccentricity(self):
        """
           Test orbit type generated from eccentricity
        """

        assert (Orbit(bodies.earth, e=0.1).type() == OrbitType.Elliptic)
        assert (Orbit(bodies.earth, e=1).type() == OrbitType.Parabolic)
        assert (Orbit(bodies.earth, e=1.1).type() == OrbitType.Hyperbolic)

    def test_from_ra(self):
        """
            Test apoapsis implies eccentric orbit
        """

        assert (Orbit(bodies.earth, ra=20000 * units.km).type() == OrbitType.Elliptic)

    def test_from_a(self):
        """
            Test orbit type detection from semimajor axis (based on sign)
        """
        assert (Orbit(bodies.earth, a=18000*units.km).type() == OrbitType.Elliptic)
        assert (Orbit(bodies.earth, a=np.inf*units.km).type() == OrbitType.Parabolic)
        assert (Orbit(bodies.earth, a=-18000*units.km).type() == OrbitType.Hyperbolic)

