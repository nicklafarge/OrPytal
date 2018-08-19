from abc import ABCMeta, abstractmethod
from utils import copydoc
import logging
import seaborn as sns
import numpy as np
import frames

class PlotUtilsBase(metaclass=ABCMeta):
    def new(self, **kwargs):
        """
        Return a new instance of the plotting engine
        """
        obj = self.__new__(self.__class__)
        obj.__init__(**kwargs)
        return obj

    def plot_traj_list(self, traj_list, **kwargs):
        """
        Plot a list of trajectories (in three-dimensions)
        :param traj_list: list of cr3bp.Trajectory
        :param kwargs: arguments to pass the the plotting function
        :return:
        """
        lines = []

        for i, traj in enumerate(traj_list):
            if (i + 1) % 10 == 0:
                logging.info("Traj {} / {}".format(i + 1, len(traj_list)))

            if 'ls' not in kwargs:
                kwargs['ls'] = 'solid'

            traj_plot = self.plot_traj(traj, **kwargs)
            lines.append(traj_plot)

        return lines

    def update_plot(self):
        """
        Update the plot that is already being shown.
        :return:
        """
        pass

    def _to_list(self, val):
        """
        Helper method that will convert a value to a list (many plotting function require [x], [y], [z]
        """
        if not hasattr(val, '__len__'):
            val = [val]
        return val

    @abstractmethod
    def init_cr3bp_plot(self):
        """
        Initialize the 3d plot. This is specific to the plotting engine.
        """
        return NotImplemented

    @abstractmethod
    def default_axes(self, *args, **kwargs):
        """
        Set the default axes for the plot. If a state is supplied, is is used to check if plot is dimensional

        Default: units are nondimensional, indecated by [nd]

        :param args: optional arguments for axes

        :param state: cr3bp.State to check for dimensionality
        """
        return NotImplemented

    @abstractmethod
    def plot_primary(self, system, **kwargs):
        """
        Plot P1 in a three body system.

        The system MUST define system.primary_body.radius.

        :param system: cr3bp.System that contains body data
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented

    @abstractmethod
    def plot_traj(self, traj, frame=frames.InertialFrame, **kwargs):
        """
        Plot a given trajectory trajectory.

        :param traj: trajectory.Trajectory to plot
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented

    @abstractmethod
    def plot_state(self, state, frame=frames.InertialFrame, **kwargs):
        """
        Plot an instance of state.State

        :param state: instance of trajectory.State
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented

    @abstractmethod
    def show(self):
        """
        Show the plot
        """
        return NotImplemented

    @abstractmethod
    def xlabel(self, xlabel, **kwargs):
        return NotImplemented

    @abstractmethod
    def ylabel(self, ylabel, **kwargs):
        return NotImplemented

    @abstractmethod
    def title(self, title, **kwargs):
        return NotImplemented


class PlotUtils2D(PlotUtilsBase, metaclass=ABCMeta):
    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, body, **kwargs):
        return self.plot_circle(0, 0, body.radius, label=body.name, **kwargs)

    @copydoc(PlotUtilsBase.default_axes)
    def default_axes(self, var1='x', var2='y', latex=False):

        unit1 = 'km'
        unit2 = 'km'

        if latex:
            fmt_str = '${} \; [{}]$'
        else:
            fmt_str = '{} [{}]'

        x_label = fmt_str.format(var1, unit1)
        y_label = fmt_str.format(var2, unit2)

        self.xlabel(x_label)
        self.ylabel(y_label)

    @copydoc(PlotUtilsBase.plot_traj)
    def plot_traj(self, traj, frame=frames.InertialFrame, **kwargs):
        return self.plot_traj_projection(traj, **kwargs)

    @copydoc(PlotUtilsBase.plot_state)
    def plot_state(self, state, frame=frames.InertialFrame, **kwargs):
        self.plot_projection(state, frame=frames.InertialFrame, **kwargs)

    @copydoc(PlotUtilsBase.plot_traj_list)
    def plot_traj_list_projection(self, traj_list, **kwargs):
        """
        Shadows 'plot_traj_list'
        """
        return self.plot_traj_list(traj_list, **kwargs)

    def plot_traj_projection(self, traj, *args, **kwargs):
        """
        Plot a given trajectory. Each trajectory will be projected on given axes (default is x-y)

        :param traj: Instance of cr3bp.Trajectory
        :param kwargs: arguments to pass the the plotting function
        """

        if 'ls' not in kwargs:
            kwargs['ls'] = 'solid'

        if 'label' not in kwargs and traj.name:
            kwargs['label'] = traj.name

        return self.plot_projection(traj.states, *args, **kwargs)

    def plot_projection(self, state_list, frame=frames.InertialFrame, var1='x', var2='y', poincare_map=False, **kwargs):
        """
        Plot a projection of a list of states.

        Each state will be projected on given axes (default is x-y)

        :param state_list: List of cr3bp.State objects
        :param var1: x-axis projection variable (default = 'x')
        :param var2: y-axis projection variable (default = 'y')
        :param poincare_map: true if this projection is a poincare map
        :param kwargs: arguments to pass the the plotting function
        """
        l_vars = ['x', 'y', 'z']
        ls_vars = ['dx', 'dy', 'dz']

        valid_vars = l_vars + ls_vars

        if var1 not in valid_vars:
            raise ValueError('%s not in (%s)' % (var1, ','.join(valid_vars)))

        if var2 not in valid_vars:
            raise ValueError('%s not in (%s)' % (var2, ','.join(valid_vars)))

        if not poincare_map and 'ls' not in kwargs:
            kwargs['ls'] = 'solid'

        if 'label' not in kwargs:
            kwargs['label'] = '_nolegend_'

        if not hasattr(state_list, '__len__'):
            state_list = [state_list]

        var1_axis = [getattr(st, var1) for st in state_list]
        var2_axis = [getattr(st, var2) for st in state_list]

        return self.plot(var1_axis, var2_axis, **kwargs)

    @abstractmethod
    def plot_circle(self, x, y, radius, **kwargs):
        """
        Plot a circle at (x,y) with a given radius

        :param x: x-axis coordinate that defines the center of the circle
        :param y: y-axis coordinate that defines the center of the circle
        :param radius: radius of the circle (be careful to match the dimensions of the plot)
        :param kwargs: arguments to pass the the plotting function
        :return:
        """
        return NotImplemented

    @abstractmethod
    def plot(self, x_vals, y_vals, **kwargs):
        """
        Simply plot a list of points defined by x_vals and y_vals

        :param x_vals: list of coordinates along the x-axis
        :param y_vals: list of coordinates along the y-axis
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented


class PlotUtils3D(PlotUtilsBase, metaclass=ABCMeta):
    @copydoc(PlotUtilsBase.plot_traj)
    def plot_traj(self, traj, frame=frames.InertialFrame, **kwargs):
        if 'color' in traj.metadata and 'c' not in kwargs:
            kwargs['c'] = traj.metadata['color']

        if traj.name and 'label' not in kwargs:
            kwargs['label'] = traj.name

        return self.plot3(
            [st.x for st in traj],
            [st.y for st in traj],
            [st.z for st in traj],
            **kwargs
        )

    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, body, **kwargs):
        return self.plot_sphere(0, 0, 0, body.radius, label=body.name)

    @copydoc(PlotUtilsBase.default_axes)
    def default_axes(self, *_, **kwargs):
        unit = 'km'
        fmt_str = '$\hat{%s}$ [%s]'

        self.xlabel(fmt_str % ('x', unit), latex=True)
        self.ylabel(fmt_str % ('y', unit), latex=True)
        self.zlabel(fmt_str % ('z', unit), latex=True)

    @copydoc(PlotUtilsBase.plot_state)
    def plot_state(self, state, frame=frames.InertialFrame, **kwargs):
        return self.scatter3_states([state], **kwargs)

    def scatter3_states(self, states, frame=frames.InertialFrame, **kwargs):
        """
        Scatter plot a list of cr3bp.State
        :param states: list of cr3bp.State
        :param kwargs: arguments to pass the the plotting function
        :return:
        """
        return self.scatter3(
            [st.x for st in states],
            [st.y for st in states],
            [st.z for st in states],
            **kwargs
        )

    @abstractmethod
    def plot3(self, x_vals, y_vals, z_vals, **kwargs):
        """
        Plot a list of x,y,z values in 3D, and connect the points (as in a trajectory)
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented

    @abstractmethod
    def scatter3(self, x_vals, y_vals, z_vals, **kwargs):
        """
        Scatter plot a list of x,y,z values in 3D
        :param kwargs: arguments to pass the the plotting function
        """
        return NotImplemented

    @abstractmethod
    def plot_sphere(self, x, y, z, radius, **kwargs):
        """
        Plot a sphere centered at (x,y,z) with a given radius
        """
        return NotImplemented
