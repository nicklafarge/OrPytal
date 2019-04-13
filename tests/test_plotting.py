########### Standard ###########
import unittest
import logging

########### Local ###########
from orpytal import Orbit, KeplarianState, plotting, frames
from orpytal.common import units
from orpytal.planet_constants import earth

########### External ###########
import numpy as np

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

# Create an orbit
orbit = Orbit(earth)
state = KeplarianState(orbit)

orbit.e = 0.1
orbit.p = 18400000 * units.m
orbit.raan = 10 * units.deg
orbit.i = 75 * units.deg

state.r = 18500 * units.km
state.arg_latitude = 55 * units.deg
state.ascending = True

traj = orbit.propagate_full_orbit()

def test_plotly_3d():
    # plotting.plot_orbit(orbit)
    plotting.plot_orbit(orbit, frame=frames.PerifocalFrame)


if __name__ == '__main__':
    test_plotly_3d()