########### Standard ###########
import logging
import itertools

########### Local ###########
import orpytal as op

# Compute an orbit given radius of apoapsis and semilatus rectum
example_orbit = op.Orbit(op.bodies.earth)

example_orbit.period = 7.3 * op.units.hours
example_orbit.e = 0.5
example_orbit.raan = 120 * op.units.deg
example_orbit.arg_periapsis = 70 * op.units.deg
example_orbit.i = 45 * op.units.deg

op.plotting.plot_orbit(example_orbit, frame=op.frames.InertialFrame)