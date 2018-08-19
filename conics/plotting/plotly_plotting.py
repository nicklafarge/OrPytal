########### Standard ###########
import logging
import itertools
import json
import re

########### Local ##########
from common import units, copydoc
from utils import copydoc
from plotting_base import PlotUtils3D, PlotUtils2D, PlotUtilsBase

########### External ###########
from dash import Dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output
from flask import Flask, render_template
import numpy as np
import plotly
import plotly.graph_objs as go
import seaborn as sns

plotly.tools.set_credentials_file(username='nlafarge', api_key='FW8T5gLFcKYHcT1fICQa', stream_ids=['qpoh56k49e'])
app = Dash(__name__)


# app = Flask(__name__)
# app.debug = True


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
        print(body_name)
        if body_name.lower() == 'moon':
            return self.moon_color
        elif body_name.lower() == 'earth':
            return self.earth_color
        elif body_name.lower() == 'sun':
            return self.sun_color
        else:
            return 'black'

    def plotly_title(self, title, **kwargs):
        self._layout.title = title

    def plotly_xlabel(self, xlabel, **kwargs):
        self._layout.xaxis.title = xlabel

    def plotly_ylabel(self, ylabel, **kwargs):
        self._layout.xaxis.title = ylabel

    def zlabel(self, zlabel):
        self._layout.xaxis.title = zlabel

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
        # self.app.run_server(debug=True)

    def xlim(self, lims):
        self._layout.xaxis.range = lims

    def ylim(self, lims):
        self._layout.yaxis.range = lims

    def zlim(self, lims):
        self._layout.zaxis.range = lims

    def plotly_update_plot(self):
        pass

        # @app.callback(Output('plot', 'children'),
        #               Input('plot', 'children'))
        # def select_points(self, value):
        #     print('--------- SELECT POINTS ---------')
        #     pass


class PlotlyPlotUtils3D(PlotUtils3D, PlotlyPlotUtils):
    def __init__(self, **kwargs):
        self.init_plotly(**kwargs)

        camera = dict(
            up=dict(x=1, y=0, z=0),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=0, y=0, z=1)
        )

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
            aspectratio={'x': 1, 'y': 1, 'z': 1},
            camera=camera
        )

    @copydoc(PlotUtilsBase.show)
    def show(self, **kwargs):
        super().show_plotly(**kwargs)

    @copydoc(PlotUtilsBase.init_cr3bp_plot)
    def init_cr3bp_plot(self, fig_number=None):
        self.__init__()

    @copydoc(PlotUtilsBase.plot_li)
    def plot_li(self, system, li, **kwargs):
        kwargs = self.format_args(**kwargs)
        li_pos = getattr(system, li.lower()).pos()

        size = 2 if 'size' not in kwargs else kwargs.pop('size')
        symbol = 'x' if 'symbol' not in kwargs else kwargs.pop('symbol')
        color = 'black' if 'color' not in kwargs else kwargs.pop('color')
        marker_dict = {
            'symbol': symbol,
            'size': size,
            'opacity': 0.7,
            'color': color
        }

        if 'name' not in kwargs:
            kwargs['name'] = li.upper()

        trace = go.Scatter3d(
            x=self._to_list(li_pos[0]), y=self._to_list(li_pos[1]), z=self._to_list(li_pos[2]),
            mode="markers",
            marker=marker_dict,
            **kwargs
        )
        self._data.append(trace)

    @copydoc(PlotUtilsBase.plot_primary)
    def plot_primary(self, system, **kwargs):
        kwargs = self.format_args(**kwargs)
        if 'name' not in kwargs:
            kwargs['name'] = system.primary_body.name
        super().plot_p1(system, **kwargs)

    @copydoc(PlotUtilsBase.plot_p2)
    def plot_p2(self, system, **kwargs):
        kwargs = self.format_args(**kwargs)
        if 'name' not in kwargs:
            kwargs['name'] = system.secondary_body.name
        super().plot_p2(system, **kwargs)

    @copydoc(PlotUtils3D.plot_sphere)
    def plot_sphere(self, x, y, z, radius, num=20, **kwargs):
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
            'width': 5 if 'width' not in kwargs else kwargs.pop('width'),
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
        self.plotly_xlabel(xlabel)

    @copydoc(PlotUtils2D.ylabel)
    def ylabel(self, ylabel, **kwargs):
        self.plotly_ylabel(ylabel)

    @copydoc(PlotUtilsBase.update_plot)
    def update_plot(self):
        self.plotly_update_plot()


class PlotlyPlotUtils2D(PlotUtils2D, PlotlyPlotUtils):
    def __init__(self):
        self.init_plotly()

    @copydoc(PlotUtilsBase.show)
    def show(self, **kwargs):
        self.show_plotly(**kwargs)

    @copydoc(PlotUtilsBase.init_cr3bp_plot)
    def init_cr3bp_plot(self):
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

    def default_axes(self, var1='x', var2='y', **kwargs):
        xlabel, ylabel = self.get_default_axes(var1, var2, latex=False, **kwargs)
        self._layout.xaxis.title = xlabel
        self._layout.yaxis.title = ylabel

    @copydoc(PlotUtils2D.plot_circle)
    def plot_circle(self, x, y, radius, num=500, **kwargs):
        kwargs = self.format_args(**kwargs)

        u1 = np.linspace(0, 2 * np.pi, num)

        xx = x + radius * np.cos(u1)
        yy = y + radius * np.sin(u1)

        color = 'black'
        if 'color' in kwargs:
            color = kwargs.pop('color')

        line = line = dict(color=color, width=5, dash='dash')
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
    def plot_primary(self, system, **kwargs):
        kwargs = self.format_args(**kwargs)
        if 'name' not in kwargs:
            kwargs['name'] = system.primary_body.name
        super().plot_p1(system, **kwargs)

    @copydoc(PlotUtilsBase.plot_p2)
    def plot_p2(self, system, **kwargs):
        kwargs = self.format_args(**kwargs)
        if 'name' not in kwargs:
            kwargs['name'] = system.secondary_body.name
        super().plot_p2(system, **kwargs)

    @copydoc(PlotUtilsBase.plot_li)
    def plot_li(self, system, li, **kwargs):
        kwargs = self.format_args(**kwargs)
        li_pos = getattr(system, li.lower()).pos()

        size = 2 if 'size' not in kwargs else kwargs.pop('size')
        symbol = 'x' if 'symbol' not in kwargs else kwargs.pop('symbol')
        color = 'black' if 'color' not in kwargs else kwargs.pop('color')
        marker_dict = {
            'symbol': symbol,
            'size': size,
            'opacity': 0.7,
            'color': color
        }

        if 'name' not in kwargs:
            kwargs['name'] = li.upper()

        trace = go.Scatter(
            x=self._to_list(li_pos[0]), y=self._to_list(li_pos[1]),
            mode="markers",
            marker=marker_dict,
            **kwargs
        )
        self._data.append(trace)

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

    @copydoc(PlotUtils2D.plot_poincare)
    def plot_poincare(self, x_vals, y_vals, **kwargs):
        kwargs = self.format_args(**kwargs)

        size = 2 if 'size' not in kwargs else kwargs.pop('size')
        symbol = 'circle' if 'symbol' not in kwargs else kwargs.pop('symbol')
        marker_dict = {
            'symbol': symbol,
            'size': size,
            'opacity': 1
        }
        if 'color' in kwargs:
            marker_dict['color'] = kwargs.pop('color')

        trace = go.Scatter(
            x=self._to_list(x_vals), y=self._to_list(y_vals),
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
        self.plotly_xlabel(xlabel)

    @copydoc(PlotUtils2D.ylabel)
    def ylabel(self, ylabel, **kwargs):
        self.plotly_ylabel(ylabel)

    @copydoc(PlotUtilsBase.update_plot)
    def update_plot(self):
        self.plotly_update_plot()
