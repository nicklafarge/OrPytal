########### Standard ###########
import logging
import unittest

########### Local ###########
from conics import Orbit, KeplarianState, plotting, frames, bodies
from common import units, Q_

########### External ###########
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=4)

orbit = Orbit(bodies.earth, name='Orbit1')
state = KeplarianState(orbit, name='State1')

orbit2 = Orbit(bodies.earth, name='Orbit1 Check')
state2 = KeplarianState(orbit2, name='State1 Check')

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

if __name__ == '__main__':
    # state2.position = frames.Vector(orbit2, state, np.array([106780, -152498, 0])*units.km, frames.InertialFrame)
    # state2.velocity = frames.Vector(orbit2, state, np.array([0.4102, 0.5637, 0])*units('km/s'), frames.InertialFrame)

    orbit.e = 0.8
    orbit.p = 37800 * units.km
    orbit.inclination = 2 * units.deg
    orbit.ascending_node = 75 * units.deg
    orbit.arg_periapsis = 30 * units.deg
    state.ta = -175 * units.deg

    print('------------ Orbit 1 ------------')
    print(orbit)

    print('------------ State 1 ------------')
    print(state)

    state2.position = frames.Vector.from_vector(orbit2, state2, state.position.inertial())
    state2.velocity = frames.Vector.from_vector(orbit2, state2, state.velocity.inertial())
    print('------------ Orbit 2 (b) ------------')
    print(orbit2)

    plt.figure(1)
    frame = 'inertial'
    plotting.plot_orbit(orbit, frame=frame)
    plotting.plot_state(state, frame=frame)
    plt.show(block=False)

    state.position.inertial()
    # plt.figure(2)
    # plotting.animate_orbit_fixed(orbit)

    pass
