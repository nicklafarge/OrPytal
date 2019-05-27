########### Standard ###########
import logging
from abc import ABCMeta, abstractmethod

########### Local ###########
from orpytal.errors import ParameterUnavailableError
from orpytal import frames, Orbit, Trajectory
from orpytal import planet_constants, units
from orpytal.utils.utils import copydoc

########### External ###########
import numpy as np


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
        Plot a list of trajectories
        :param traj_list: list of trajectory.Trajectory objects
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
        """
        pass

    def _to_list(self, val):
        """
        Helper method that will convert a value to a list (many plotting function require [x], [y], [z]
        """
        if not hasattr(val, '__len__'):
            val = [val]
        return val

    def init_plot(self, frame=frames.InertialFrame,):
        """
        Initialize the plot.

        :param frame: the frame the plot will be in
        """
        self.frame = frame

    def plot_traj(self, traj, **kwargs):
        """
        Plot a given trajectory.

        :param traj: trajectory.Trajectory to plot
        :param kwargs: arguments to pass the the plotting function
        """

        traj.in_frame(self.frame)

    def plot_state(self, state, **kwargs):
        """
        Plot an orbital State

        :param state: instance of state.State
        :param kwargs: arguments to pass the the plotting function
        """
        state.to(self.frame, state)

    @abstractmethod
    def default_axes(self, *args, **kwargs):
        """
        Set the default axes for the plot. If a state is supplied, is is used to check if plot is dimensional
        """
        return NotImplemented


    @abstractmethod
    def plot_primary(self, primary, **kwargs):
        """
        Plot The primary body

        :param primary: some object that contains primary body data
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
    def plot_primary(self, primary, **kwargs):

        if isinstance(primary, Orbit) or isinstance(primary, Trajectory):
            radius = primary.central_body.radius
        elif isinstance(primary, planet_constants.CentralBody):
            radius = primary.radius
        else:
            radius = primary

        primary_location = [0, 0]

        if isinstance(radius, units.Quantity):
            radius = radius.m

        return self.plot_circle(*primary_location, radius)


    @copydoc(PlotUtilsBase.default_axes)
    def default_axes(self,frame=frames.InertialFrame):

        unit = 'km'
        fmt_str = '%s[%s]'
        labels = ['e', 'p']

        self.xlabel(fmt_str % (labels[0], unit), latex=True)
        self.ylabel(fmt_str % (labels[1], unit), latex=True)

    @copydoc(PlotUtilsBase.plot_traj)
    def plot_traj(self, traj, **kwargs):
        super().plot_traj(traj, **kwargs)

        if 'ls' not in kwargs:
            kwargs['ls'] = 'solid'

        if 'label' not in kwargs and traj.name:
            kwargs['label'] = traj.name

        return self.plot(traj.x_vals, traj.y_vals)

    @copydoc(PlotUtilsBase.plot_state)
    def plot_state(self, state, **kwargs):
        state.to(self.frame, state)
        self.plot_projection(state, **kwargs)

    @copydoc(PlotUtilsBase.plot_traj_list)
    def plot_traj_list_projection(self, traj_list, **kwargs):
        """
        Shadows 'plot_traj_list'
        """
        return self.plot_traj_list(traj_list, **kwargs)


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
    def plot_traj(self, traj, **kwargs):
        super().plot_traj(traj, **kwargs)

        if 'color' in traj.metadata and 'c' not in kwargs:
            kwargs['c'] = traj.metadata['color']

        if traj.name and 'label' not in kwargs:
            kwargs['label'] = traj.name


        return self.plot3(traj.x_vals, traj.y_vals, traj.z_vals, **kwargs)

    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, primary, **kwargs):

        if isinstance(primary, Orbit) or isinstance(primary, Trajectory):
            radius = primary.central_body.radius
        elif isinstance(primary, planet_constants.CentralBody):
            radius = primary.radius
        else:
            radius = primary

        primary_location = [0, 0, 0]

        if isinstance(radius, units.Quantity):
            radius = radius.m

        return self.plot_sphere(*primary_location, radius, **kwargs)

    @copydoc(PlotUtilsBase.default_axes)
    def default_axes(self, frame=frames.InertialFrame):
        unit = 'km'
        fmt_str = '{}[{}]'

        if frame == frames.InertialFrame:
            labels = ['x', 'y', 'z']
        elif frame == frames.PerifocalFrame:
            labels = ['e', 'p', 'h']
        else:
            raise ValueError(f"Can't set axes for {frame.name}")



        self.xlabel(fmt_str.format(labels[0], unit))
        self.ylabel(fmt_str.format(labels[1], unit))
        self.zlabel(fmt_str.format(labels[2], unit))

    @copydoc(PlotUtilsBase.plot_state)
    def plot_state(self, state, **kwargs):
        return self.scatter3_states([state], **kwargs)

    def scatter3_states(self, states, **kwargs):
        """
        Scatter plot a list of States
        :param states: list of state.State
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
