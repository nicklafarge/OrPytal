########### Standard ###########

########### Local ###########
from conics import Orbit, KeplarianState

########### External ###########
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def plot_orbit_fixed(orbit, step=0.01, **kwargs):
    ta_range = np.arange(0, 2*np.pi, step)

    x_list = []
    y_list = []

    for ta in ta_range:
        st = KeplarianState(orbit)
        st.ta = ta
        orbit_fixed_pos = st.position.orbit_fixed()
        x_list.append(orbit_fixed_pos[0])
        y_list.append(orbit_fixed_pos[1])

    plt.axis('equal')
    plt.plot(x_list, y_list, **kwargs)

    plot_primary(orbit)
    plt.xlabel("{} \; [{}]".format('$\hat{e}$', 'km'))
    plt.ylabel("{} \; [{}]".format('$\hat{p}$', 'km'))


def plot_primary(orbit):
    fig = plt.gcf()
    ax = fig.gca()

    circle_primary = plt.Circle((0, 0), orbit.central_body.radius.m,
                                color='k',
                                fill=False,
                                lw=2.5)
    ax.add_artist(circle_primary)
