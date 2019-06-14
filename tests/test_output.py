########### Standard ###########
import itertools
import logging
import unittest

########### Local ###########
from orpytal import Orbit, KeplarianState, frames, bodies, OrbitType
from orpytal.common import units
from orpytal.planet_constants import CentralBody

########### External ###########
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit as PoliastroOrbit
from astropy import units as u

a = 51000 * units.km
i = 1.85 * units.deg
raan = 49.562 * units.deg
argp = 286.537 * units.deg

elliptic_orbit = Orbit(bodies.earth, a=a, e=0.1, raan=raan, arg_periapsis=argp, i=i, name='TestOrbit')
hyperbolic_orbit = Orbit(bodies.earth, a=-a, e=1.1, raan=raan, arg_periapsis=argp, i=i)


class TestOrbitCreation(unittest.TestCase):
    """
        Unit tests for orbit creation (circular/elliptic)
    """

    def test_elliptic(self):
        orbit_output = str(elliptic_orbit)
        assert ('Orbit Name: TestOrbit' in orbit_output)
        assert ('Rad. Apoapsis' in orbit_output)
        assert ('Orbital Period' in orbit_output)

    def test_hyperbolic(self):
        orbit_output = str(hyperbolic_orbit)

        # Hyperbolic Quantities
        assert ('V Infinity' in orbit_output)
        assert ('TA Infinity' in orbit_output)
        assert ('Flyby Angle' in orbit_output)

    def test_state_full(self):
        st = elliptic_orbit.get_state(ta=60 * units.deg, name='State1')
        state_output = str(st)

        assert('State Name: State1' in state_output)
        assert('Inertial Frame' in state_output)
        assert(' x: ' in state_output)
        assert(' y: ' in state_output)
        assert(' z: ' in state_output)

    def test_state_perifocal(self):

        elliptic_orbit_perifocal = Orbit(bodies.earth, a=a, e=0.1)
        st = elliptic_orbit_perifocal.get_state(ta=60 * units.deg)
        state_output = str(st)

        # assert('State Name: State1' in state_output)
        assert('Perifocal Frame' in state_output)
        assert(' e: ' in state_output)
        assert(' p: ' in state_output)
        assert(' h: ' in state_output)

    def test_state_perifocal_no_orbit_info(self):

        elliptic_orbit_perifocal = Orbit(bodies.earth)
        st = elliptic_orbit_perifocal.get_state(r=26000*units.km, ta=60 * units.deg)
        state_output = str(st)

        assert('Perifocal Frame' in state_output)
        assert(' e: ' in state_output)
        assert(' p: ' in state_output)
        assert(' h: ' in state_output)

    def test_state_rotating(self):

        elliptic_orbit_perifocal = Orbit(bodies.earth)
        st = elliptic_orbit_perifocal.get_state(r=26000*units.km)
        state_output = str(st)

        assert('Rotating Frame' in state_output)
        assert(' r: ' in state_output)
        assert(' theta: ' in state_output)
        assert(' h: ' in state_output)

    def test_state_no_information(self):

        elliptic_orbit_perifocal = Orbit(bodies.earth)
        st = elliptic_orbit_perifocal.get_state()
        state_output = str(st)

        assert('Central Body: Earth' in state_output)
        assert('Inertial Frame' in state_output)