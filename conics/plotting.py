########### Standard ###########
import logging

########### Local ###########
from conics import Orbit, KeplarianState, state
from common import units, Q_

########### External ###########
import matplotlib.animation as animation
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


def plot_orbit_fixed(orbit, **kwargs):
    return plot_orbit(orbit, frame='orbit_fixed', **kwargs)


def plot_orbit_inertial(orbit, **kwargs):
    return plot_orbit(orbit, frame='inertial', **kwargs)


def plot_orbit_inertial_3d(orbit, **kwargs):
    return plot_orbit(orbit, frame='inertial', planar=False, **kwargs)


def plot_orbit(orbit, frame='orbit_fixed', planar=True, **kwargs):
    logging.info('Propagating and plotting {}'.format(orbit.name))

    start_st = KeplarianState(orbit)
    start_st.ta = 0 * units.rad
    traj = orbit.propagate_full_orbit(start_st)
    x_list, y_list, z_list = traj.inertial()

    if 'label' not in kwargs:
        kwargs['label'] = orbit.name

    logging.info('Propagation Complete')

    if planar:
        plt.axis('equal')
        plt.plot(x_list, z_list, ls='dashed', **kwargs)
        plot_primary(orbit)
        plt.title(orbit.name)
        if frame == 'orbit_fixed':
            xaxis = 'e'
            yaxis = 'p'
        elif frame == 'inertial':
            xaxis = 'x'
            yaxis = 'y'
        else:
            return

        plt.xlabel("{} [{}]".format('$\hat{%s}$' % xaxis, 'km'))
        plt.ylabel("{} [{}]".format('$\hat{%s}$' % yaxis, 'km'))
        plt.legend()
    else:
        import matlab_plotting as mplt
        mplt.init_cr3bp_3d_plot(1)
        mplt.title(orbit.name)
        mplt.default_3d_axes()
        mplt.plot_3d_traj(traj, label=orbit.name)
        mplt.show()


def plot_state(state):
    orbit_fixed_pos = state.position.orbit_fixed()
    plt.plot(orbit_fixed_pos[0], orbit_fixed_pos[1],
             marker='o',
             c='r')


def plot_primary(orbit):
    fig = plt.gcf()
    ax = fig.gca()

    circle_primary = plt.Circle((0, 0), orbit.central_body.radius.m,
                                color='k',
                                fill=True,
                                lw=2.5)
    ax.add_artist(circle_primary)


def animate_orbit_fixed(orbit, step=0.01, **kwargs):
    ta_range = np.arange(0, 2 * np.pi, step)

    x_list = []
    y_list = []

    logging.info('Propagating... ')
    traj = []
    for ta in ta_range:
        st = KeplarianState(orbit)
        st._ta.value = ta * units.radians
        st._r.set(st, orbit)
        st._position.set(st, orbit)
        orbit_fixed_pos = st.position.orbit_fixed()
        x_list.append(orbit_fixed_pos[0])
        y_list.append(orbit_fixed_pos[1])
        traj.append(st)

    logging.info('Propagation Complete')
    traj_animator = TrajectoryAnimator(orbit, traj)
    traj_animator.animate()


class TrajectoryAnimator(object):
    def __init__(self,
                 orbit,
                 traj,
                 ta_fmt_str="ta = %0.4f [rad]"):

        self.traj = traj
        self.ta_fmt_str = ta_fmt_str

        plot_primary(orbit)
        plt.axis("equal")
        plt.title(orbit.name)
        plt.xlabel("{} [{}]".format('$\hat{e}$', 'km'))
        plt.ylabel("{} [{}]".format('$\hat{p}$', 'km'))
        ax = plt.gcf().gca()

        self.ta_text = ax.text(.05, 0.05, '',
                               fontsize=15,
                               horizontalalignment='left',
                               verticalalignment='bottom',
                               transform=plt.gca().transAxes)

        self.t_loc = ax.plot([], [], 'bo')[0]

        # For the axes....
        self.orbit_fixed_pos = [st.position.orbit_fixed() for st in self.traj]
        plt.plot([st[0] for st in self.orbit_fixed_pos],
                 [st[1] for st in self.orbit_fixed_pos],
                 c='0.85', ls='dashed')

    def _init(self):
        return self._animate(0)

    def _animate(self, i):
        st = self.traj[i]

        self.t_loc.set_data(self.orbit_fixed_pos[i][0], self.orbit_fixed_pos[i][1])

        if self.ta_fmt_str:
            self.ta_text.set_text(self.ta_fmt_str % self.traj[i].ta.m)

        return self.t_loc, self.ta_text

    def animate(self, save_as_gif_name=None, interval=20, fps=1):
        frames = np.arange(0, len(self.traj))
        ani = animation.FuncAnimation(plt.gcf(), self._animate, frames,
                                      interval=interval, blit=True)
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

        if save_as_gif_name:
            print('Saving gif...')
            ani.save(save_as_gif_name, writer='imagemagick', fps=fps)
            print('Saving complete!')
        else:
            print('Animating...')
            plt.show()
