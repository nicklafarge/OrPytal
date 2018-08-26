########### Standard ###########
import logging

########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting

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
    # state.position = frames.Vector(orbit2, state, np.array([106780, -152498, 0])*units.km, frames.InertialFrame)
    # state.velocity = frames.Vector(orbit2, state, np.array([0.4102, 0.5637, 0])*units('km/s'), frames.InertialFrame)

    orbit.e = 0.8
    orbit.p = 37800 * units.km
    orbit.inclination = 45 * units.deg
    orbit.ascending_node = 75 * units.deg
    orbit.arg_periapsis = 30 * units.deg
    # state.ta = -175 * units.deg

    # state.r = 37800 * units.km
    print('------------ Orbit 1 ------------')
    print(orbit)

    print('------------ State 1 ------------')
    print(state)

    orbit2.e = orbit.e
    orbit2.p = orbit.p
    orbit2.inclination = orbit.inclination
    orbit2.arg_periapsis = orbit.arg_periapsis
    orbit2.ascending_node = 120 * units.deg
    plt.figure(1)
    frame = 'inertial'
    plotting.plot_orbit(orbit, frame=frame, planar=False)
    plotting.plot_orbit(orbit2, frame=frame, planar=False)
    plotting.plot_primary(orbit)
    # plotting.plot_state(state, frame=frame)
    plt.show(block=False)
    #
    # state.position.inertial()
    # plt.figure(2)
    # plotting.animate_orbit_fixed(orbit)

    pass
