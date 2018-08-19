########### Standard ###########
import logging

########### Local ###########
from cr3bp_plotting.plotting_base import PlotUtils3D
from utils import utils, cr3bp_utils

########### External ###########
import matlab.engine
import numpy as np
import seaborn as sns

logging.info("Starting Matlab Engine....")
eng = matlab.engine.start_matlab()
logging.info("Matlab Engine Started!")

class MatlabPlotUtils3D(PlotUtils3D):
    def __init__(self):
        self.eng = eng

    def eval(self, cmd, nargout=0):
        self.eng.eval(cmd, nargout=nargout)

    def figure(self, fig_number=None):
        if fig_number:
            self.eng.figure(float(fig_number))
        else:
            self.eng.figure()

    def init_cr3bp_plot(self, fig_number=None):
        self.figure(fig_number)

        self.eval("hold on;")
        self.eval("axis equal;")
        self.default_axes()

    def plot_sphere(self, x, y, z, radius, **kwargs):
        self.eng.workspace['radius'] = radius
        self.eng.workspace['x'] = float(x)
        self.eng.workspace['y'] = float(y)
        self.eng.workspace['z'] = float(z)
        self.eval("[sx sy sz] = sphere;")
        self.eval("surf(sx*radius + x, sy*radius + y, sz*radius + z, 'FaceColor', [0 0 0]);")

    def legend(self, *args):
        if not args:
            self.eng.legend()
        else:
            self.eng.legend(*args)

    def _prepare_plot(self, x_vals, y_vals, z_vals, **kwargs):
        x_vals = self._to_list(x_vals)
        y_vals = self._to_list(y_vals)
        z_vals = self._to_list(z_vals)

        x_vals = matlab.double(x_vals)
        y_vals = matlab.double(y_vals)
        z_vals = matlab.double(z_vals)
        args = self.format_args(**kwargs)
        return x_vals, y_vals, z_vals, args

    def scatter3(self, x_vals, y_vals, z_vals, **kwargs):
        x_vals, y_vals, z_vals, args = self._prepare_plot(x_vals, y_vals, z_vals, **kwargs)

        args = self.replace_arg_key(args, 'Color', 'MarkerEdgeColor')

        return self.eng.scatter3(x_vals, y_vals, z_vals, *args, nargout=1)

    def replace_arg_key(self, args, key, new_key):
        args_list = list(args)
        for i, arg in enumerate(args_list):
            if arg == key:
                args_list[i] = new_key
        return tuple(args_list)

    def format_args(self, **kwargs):
        args = ()
        if 'c' in kwargs:
            color = kwargs['c']
            if isinstance(color, str):
                pass  # no op
            elif hasattr(color, '__len__'):
                color = matlab.double(list(color))[0]
            else:
                raise TypeError('need list or string for color')

            args = args + ("Color", color)

        if 'ls' in kwargs:
            ls = kwargs['ls']
            if ls == 'dashed':
                ls = '--'
            elif ls == 'dotted':
                ls = ':'
            args = args + ('LineStyle', ls)

        if 'label' in kwargs:
            args = args + ('DisplayName', kwargs['label'])

        if 'marker' in kwargs:
            args = args + ("Marker", kwargs['marker'])

        return args

    def plot3(self, x_vals, y_vals, z_vals, **kwargs):
        x_vals, y_vals, z_vals, args = self._prepare_plot(x_vals, y_vals, z_vals, **kwargs)

        lw = 2 if 'lw ' not in kwargs else kwargs['lw']
        args = args + ("LineWidth", lw)

        return self.eng.plot3(x_vals, y_vals, z_vals, *args, nargout=1)

    def title(self, title, latex=True):
        args = ()
        if latex:
            args = args + ("Interpreter", 'Latex')
        return self.eng.title(title, *args)

    def xlabel(self, xlabel, latex=True):
        args = ()
        if latex:
            args = args + ("Interpreter", 'Latex')

        return self.eng.xlabel(xlabel, *args)

    def ylabel(self, ylabel, latex=True):
        args = ()
        if latex:
            args = args + ("Interpreter", 'Latex')

        return self.eng.ylabel(ylabel, *args)

    def zlabel(self, zlabel, latex=True):
        args = ()
        if latex:
            args = args + ("Interpreter", 'Latex')

        return self.eng.zlabel(zlabel, *args)

    def grid(self):
        self.eng.eval("grid on;", nargout=0)

    def _axis_grid(self, axis, base, frequency, additional_ticks=10):
        axis = axis.lower()

        if axis not in ['x', 'y', 'z']:
            raise ValueError('Axis must be x, y, or z')

        limits = self.eng.eval("{}lim".format(axis))[0]

        # ax_min = utils.round_to_base(limits[0], base) - base / 2. - additional_ticks*base
        # ax_max = utils.round_to_base(limits[1], base) - base / 2. + additional_ticks*base
        ax_min = utils.round_to_base(limits[0], base) - base / 2.
        ax_max = utils.round_to_base(limits[1], base) - base / 2.

        tick_location_list = np.arange(ax_min, ax_max, step=base).tolist()
        self.eng.workspace['{}_ticks'.format(axis)] = matlab.double(tick_location_list)

        self.eval("set(gca,'{}tick',{}_ticks)".format(axis, axis))

        # self.eng.workspace['{}_label_frequency'.format(axis)] = int(frequency)

        # # X Label Frequency
        # self.eng.eval("ax = gca;", nargout=0)
        # self.eng.eval("x_labels = string(ax.XAxis.TickLabels);", nargout=0)
        # self.eng.eval("x_labels(x_label_frequency:x_label_frequency:end) = nan;", nargout=0)
        # self.eng.eval("ax.XAxis.TickLabels = x_labels;", nargout=0)

    def digraph_grid(self, base, x_frequency=2, y_frequency=2, z_frequency=2):
        self.eng.workspace['base'] = base

        self._axis_grid('x', base, x_frequency)
        self._axis_grid('y', base, y_frequency)
        self._axis_grid('z', base, z_frequency)

        # Show Grid
        self.grid()

    def adjust_tick_label_frequency(self, axis, freq=2):
        axis = axis.lower()

        if axis not in ['x', 'y', 'z']:
            raise ValueError('Axis must be x, y, or z')

        self.eng.workspace['freq'] = int(freq)

        self.eval('ax = gca;')
        # self.eval("xti = xticks;", nargout=0)
        # self.eval("ax.XAxis.TickLabels = xti;", nargout=0)

        # self.eval('labels = string({}ticks);'.format(axis))
        self.eval('labels = string(ax.{}Axis.TickLabels);'.format(axis.upper()))

        self.eval('new_labels = strings(length(labels), 1);')
        self.eval(
            'for i=1:length(labels); if mod(i,freq)==0; new_labels(i) = labels(i); else; new_labels(i)=nan; end; end;')
        self.eval('ax.{}Axis.TickLabels = new_labels;'.format(axis.upper()))
        pass

    def plot_li(self, system, li, **kwargs):
        return NotImplemented

    def show(self):
        pass
