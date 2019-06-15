########### Standard ###########
import itertools

########### Local ##########
from orpytal.plotting.plotting_base import PlotUtils3D, PlotUtils2D, PlotUtilsBase
from orpytal import Orbit, planet_constants, Trajectory, frames
from orpytal.utils.utils import copydoc

########### External ###########
from dash import Dash
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import plotly
import plotly.plotly as py
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go
import seaborn as sns

app = Dash(__name__)

BODY_COLORS = dict(
    moon='rgb(162,168,174)',
    earth="#204a87",
    sun='rgb(253,184,19)'
)


class PlotlyPlotUtils(object):
    moon_color = 'rgb(162,168,174)'
    earth_color = "#204a87"
    sun_color = 'rgb(253,184,19)'

    def init_plotly(self, filename=None):
        self.filename = filename

        self._data = []
        self._color_cycle = itertools.cycle(plotly.colors.DEFAULT_PLOTLY_COLORS)
        self._layout = go.Layout(autosize=True, height=800)

        self.app = app
        self.app.css.config.serve_locally = True
        self.app.scripts.config.serve_locally = True

        self.graph = dcc.Graph(
            id='plot',
            figure=self.figure
        )

        self.app.layout = html.Div(children=[
            dcc.Graph(
                id='plot',
                figure=self.figure
            )
        ])

    @property
    def figure(self):
        return dict(
            data=self._data,
            layout=self._layout,
        )

    def format_args(self, **kwargs):
        if 'c' in kwargs:
            c = kwargs.pop('c')
            if c == 'k':
                c = 'black'
            elif c == 'r':
                c = 'red'
            elif c == 'b':
                c = 'blue'
            elif hasattr(c, '__len__'):
                c = 'rgb({})'.format(",".join([str(i) for i in c]))
            kwargs['color'] = c

        if 'ls' in kwargs:
            ls = kwargs.pop('ls')
            if ls == 'dashed':
                ls = 'dash'
            kwargs['dash'] = ls

        if 'label' in kwargs:
            kwargs['name'] = kwargs.pop('label')

        if 'marker' in kwargs:
            s = kwargs.pop('marker')
            if s == 'o':
                s = 'circle'
            kwargs['symbol'] = s

        if 'label' in kwargs:
            kwargs['name'] = kwargs.pop('label')

        return kwargs

    def get_body_color(self, body_name):
        return BODY_COLORS.get(body_name.lower(), 'black')

    def plotly_title(self, title, **kwargs):
        self._layout.title = title

    def plotly_xlabel(self, xlabel, **kwargs):
        self._layout.xaxis.title = xlabel

    def plotly_ylabel(self, ylabel, **kwargs):
        self._layout.yaxis.title = ylabel

    def plotly_zlabel(self, zlabel, **kwargs):
        if hasattr(self._layout, 'zaxis'):
            self._layout.zaxis.title = zlabel

    def grid(self):
        pass

    def hide_legend(self):
        self._layout.showlegend = False

    def show_plotly(self, **kwargs):
        self._layout.update(**kwargs)

        plot_kwargs = {}
        if self.filename:
            plot_kwargs['filename'] = self.filename
        elif 'filename' in kwargs:
            plot_kwargs['filename'] = kwargs.pop('filename')
        plotly.offline.plot(self.figure, **plot_kwargs)

    def xlim(self, lims):
        self._layout.xaxis.range = lims

    def ylim(self, lims):
        self._layout.yaxis.range = lims

    def zlim(self, lims):
        self._layout.zaxis.range = lims

    def plotly_update_plot(self):
        pass

    def _add_primary_name(self, primary, **kwargs):
        if 'name' not in kwargs:
            if isinstance(primary, Orbit) or isinstance(primary, Trajectory):
                name = primary.central_body.name
            elif isinstance(primary, planet_constants.CentralBody):
                name = primary.name
            else:
                name = "central_body"

            if name:
                kwargs['name'] = name
        return kwargs


class PlotlyPlotUtils3D(PlotUtils3D, PlotlyPlotUtils):
    def __init__(self, frame=frames.InertialFrame):
        self.init_plotly()

        # TODO update default camera to something not cr3bp based

        if frame == frames.InertialFrame:
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.25, y=1.25, z=1.25)
            )
        elif frame == frames.PerifocalFrame:
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=-1.25, y=0, z=1.25)
            )
        else:
            raise ValueError(f"Can't plot in {self.frame.name}")

        self._layout.scene = dict(
            xaxis=dict(
                title="x [km]",
            ),
            yaxis=dict(
                title="y [km]",
            ),
            zaxis=dict(
                title="z [km]",
            ),
            aspectmode="data",
            aspectratio=dict(x=1, y=1, z=1),
            camera=camera
        )

        self.default_axes(frame=frame)

    @copydoc(PlotUtilsBase.show)
    def show(self, **kwargs):
        self.show_plotly(**kwargs)

    @copydoc(PlotUtilsBase.init_plot)
    def init_plot(self, fig_number=None, **kwargs):
        super().init_plot(**kwargs)
        self.__init__(frame=self.frame)

    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, primary, **kwargs):
        kwargs = self.format_args(**kwargs)
        kwargs = self._add_primary_name(primary, **kwargs)
        super().plot_primary(primary, **kwargs)

    @copydoc(PlotUtils3D.plot_sphere)
    def plot_sphere(self, x, y, z, radius, num=50, **kwargs):
        kwargs = self.format_args(**kwargs)
        u1 = np.linspace(0, 2 * np.pi, num)
        v1 = u1.copy()
        uu, vv = np.meshgrid(u1, v1)

        xx = x + radius * np.cos(uu) * np.sin(vv)
        yy = y + radius * np.sin(uu) * np.sin(vv)
        zz = z + radius * np.cos(vv)

        color = self.get_body_color(kwargs.pop('name')) if 'color' not in kwargs else kwargs.pop('color')

        sphere = go.Surface(
            x=xx, y=yy, z=zz,
            colorscale=[[0, color], [1, color]],
            cauto=False, cmin=1, cmax=1, showscale=False,
            **kwargs
        )

        self._data.append(sphere)

        return xx, yy, zz

    @copydoc(PlotUtils3D.plot3)
    def plot3(self, x_vals, y_vals, z_vals, **kwargs):
        kwargs = self.format_args(**kwargs)

        line_args = {
            'width': 4 if 'width' not in kwargs else kwargs.pop('width'),
            'dash': 'solid' if 'dash' not in kwargs else kwargs.pop('dash')
        }

        if 'color' in kwargs:
            line_args['color'] = kwargs.pop('color')

        trace = go.Scatter3d(
            x=x_vals, y=y_vals, z=z_vals,
            line=line_args,
            mode="lines",
            **kwargs
        )
        self._data.append(trace)

    @copydoc(PlotUtils3D.scatter3)
    def scatter3(self, x_vals, y_vals, z_vals, **kwargs):
        kwargs = self.format_args(**kwargs)

        size = 5 if 'size' not in kwargs else kwargs.pop('size')
        symbol = 'circle' if 'symbol' not in kwargs else kwargs.pop('symbol')
        marker_dict = {
            'symbol': symbol,
            'size': size,
            'opacity': 0.7
        }
        if 'color' in kwargs:
            marker_dict['color'] = kwargs.pop('color')

        trace = go.Scatter3d(
            x=self._to_list(x_vals), y=self._to_list(y_vals), z=self._to_list(z_vals),
            mode="markers",
            marker=marker_dict,
            **kwargs
        )
        self._data.append(trace)

    @copydoc(PlotUtils2D.title)
    def title(self, title, **kwargs):
        self.plotly_title(title)

    @copydoc(PlotUtils2D.xlabel)
    def xlabel(self, xlabel, **kwargs):
        self._layout.scene.xaxis.title = xlabel

    @copydoc(PlotUtils2D.xlabel)
    def ylabel(self, ylabel, **kwargs):
        self._layout.scene.yaxis.title = ylabel

    @copydoc(PlotUtils2D.ylabel)
    def zlabel(self, zlabel, **kwargs):
        self._layout.scene.zaxis.title = zlabel

    @copydoc(PlotUtilsBase.update_plot)
    def update_plot(self):
        self.plotly_update_plot()


class PlotlyPlotUtils2D(PlotUtils2D, PlotlyPlotUtils):
    def __init__(self):
        self.init_plotly()

    @copydoc(PlotUtilsBase.show)
    def show(self, **kwargs):
        self.show_plotly(**kwargs)

    @copydoc(PlotUtilsBase.init_plot)
    def init_plot(self, **kwargs):
        super().init_plot(**kwargs)

        self._layout.xaxis = dict(
            constrain="domain",
            scaleratio=1
        )

        self._layout.yaxis = dict(
            constrain="domain",
            scaleanchor="x",
            scaleratio=1
        )
        self._layout.shapes = []

        self.default_axes()

    @copydoc(PlotUtils2D.plot_circle)
    def plot_circle(self, x, y, radius, num=500, **kwargs):
        kwargs = self.format_args(**kwargs)

        u1 = np.linspace(0, 2 * np.pi, num)

        xx = x + radius * np.cos(u1)
        yy = y + radius * np.sin(u1)

        color = 'black'
        if 'color' in kwargs:
            color = kwargs.pop('color')

        line = dict(color=color, width=5, dash='dash')
        trace = go.Scatter(x=xx, y=yy, mode='markers', line=line, **kwargs)
        self._layout["shapes"] += (
            {
                'type': 'circle',
                'xref': 'x',
                'yref': 'y',
                'x0': (x - radius),
                'y0': (y - radius),
                'x1': (x + radius),
                'y1': (y + radius),
                'opacity': 1,
                'fillcolor': color,
                'line': {
                    'color': color,
                },
            },

        )
        return trace

    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, primary, **kwargs):
        kwargs = self.format_args(**kwargs)
        kwargs = self._add_primary_name(primary, **kwargs)
        super().plot_primary(primary, **kwargs)

    @copydoc(PlotUtils2D.plot)
    def plot(self, x_vals, y_vals, **kwargs):
        kwargs = self.format_args(**kwargs)

        line_args = {
            'width': 2 if 'width' not in kwargs else kwargs.pop('width'),
            'dash': 'solid' if 'dash' not in kwargs else kwargs.pop('dash')
        }

        if 'color' in kwargs:
            line_args['color'] = kwargs.pop('color')

        trace = go.Scatter(
            x=self._to_list(x_vals), y=self._to_list(y_vals),
            line=line_args,
            mode="lines",
            **kwargs
        )
        self._data.append(trace)

    @copydoc(PlotUtils2D.title)
    def title(self, title, **kwargs):
        self.plotly_title(title)

    @copydoc(PlotUtils2D.xlabel)
    def xlabel(self, xlabel, **kwargs):
        self.plotly_xlabel(xlabel)

    @copydoc(PlotUtils2D.ylabel)
    def ylabel(self, ylabel, **kwargs):
        self.plotly_ylabel(ylabel)

    @copydoc(PlotUtilsBase.update_plot)
    def update_plot(self):
        self.plotly_update_plot()
