########### Standard ###########
import logging
import itertools
########### Local ###########
from orpytal import units, Orbit, KeplarianState, bodies, plotting, get_plot_utils, integration, frames, Trajectory

import numpy as np
from matplotlib import pyplot as plt


# orbit = Orbit(bodies.earth, a=51000 * units.km, e=0.7, raan=10*units.deg, arg_periapsis=10*units.deg, inclination=45*units.deg)
orbit = Orbit(bodies.earth, a=51000 * units.km, e=0.7)
test = integration.integrate_orbit(orbit)

#
# pu = get_plot_utils("plotly", planar=False)
# pu.init_plot()
# pu.plot3(*test.y[0:3])
# pu.plot_primary(bodies.earth)
# pu.show()

r_list = [np.linalg.norm([test.y[0][i], test.y[1][i], test.y[2][i]]) for i in range(len(test.y[0]))]

initial_pos = np.array([x[0] for x in test.y[0:3]])
initial_vel = np.array([x[0] for x in test.y[3:6]])

final_state = np.array([x[-1] for x in test.y[0:3]])
error = final_state - initial_pos
print(np.linalg.norm(error))

orbit2 = Orbit(bodies.earth)
st2 = KeplarianState(orbit2)
st2.position = frames.Vector(orbit2, st2, initial_pos*units.km, frames.OrbitFixedFrame)
st2.velocity = frames.Vector(orbit2, st2, initial_vel*units("km/s"), frames.OrbitFixedFrame)
print(st2)

traj = Trajectory([KeplarianState(orbit,
                    position=frames.Vector(orbit2, st2, np.array([x[i] for x in test.y[0:3]])*units.km, frames.OrbitFixedFrame),
                    velocity=frames.Vector(orbit2, st2, np.array([x[i] for x in test.y[3:6]])*units("km/s"), frames.OrbitFixedFrame))
                   for i in range(len(test.y[0]))])
pu = get_plot_utils("plotly", planar=False)
pu.init_plot()
pu.plot3(*test.y[0:3])
pu.plot_primary(bodies.earth)
pu.show()