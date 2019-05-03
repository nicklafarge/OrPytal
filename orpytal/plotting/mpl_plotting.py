########### Standard ###########
import logging

########### Local ###########
from orpytal.plotting.plotting_base import PlotUtils2D, PlotUtils3D

########### External ###########
import matplotlib.animation as animation
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


class MplPlotUtilsBase(object):
    def plot_circle(x, y, radius, **kwargs):
        circle = plt.Circle((x, y), radius, **kwargs)
        plt.gca().add_artist(circle)

    def figure(self, **kwargs):
        return plt.figure(**kwargs)

    def plot_sphere_matplotlib(x, y, z, radius):
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = x + radius * np.outer(np.cos(u), np.sin(v))
        y = y + radius * np.outer(np.sin(u), np.sin(v))
        z = z + radius * np.outer(np.ones(np.size(u)), np.cos(v))

        plt.gca().plot_surface(x, y, z, rstride=4, cstride=4, color='r', linewidth=0, alpha=0.9)

    def plot_arrows(line, n=3, offset_pct=0.0, length=1e-4, width=1e-4, tol=1e-8, rev=False):
        ax = plt.gca()

        data = line.get_data(orig=True)
        x_vals = data[0]
        y_vals = data[1]
        c = line.get_color()

        step_size = len(x_vals) / n

        offset = int(len(x_vals) * offset_pct)

        for i in range(0, n):
            j = (i * step_size + offset) % len(x_vals)

            x = x_vals[j]
            y = y_vals[j]

            r = np.array([x, y])

            dir = 1 if not rev else -1
            dx = (x_vals[j + dir] - x)
            dy = (y_vals[j + dir] - y)

            if abs(dx) > tol or abs(dy) > tol:
                dx_vector = np.array([dx, dy])
                dx_hat = dx_vector / np.linalg.norm(dx_vector)

                # print 'arrow at (%0.5f, %0.5f) for j=%i' % (r[0], r[1], j)
                ax.arrow(r[0], r[1], dx_hat[0] * length, dx_hat[1] * length, color=c, width=width)
            else:
                print('Arrow not shown for same dx/dy')

    def custom_legend(labels, colors):
        plt.gca().legend(handles=[plt.Line2D([0], [], color=colors[i], label=labels[i]) for i in range(len(labels))])

    def mpl_title(self, title, *args, **kwargs):
        return plt.title(title, *args, **kwargs)

    def mpl_xlabel(self, title, *args, **kwargs):
        return plt.xlabel(title, *args, **kwargs)

    def mpl_ylabel(self, title, *args, **kwargs):
        return plt.ylabel(title, *args, **kwargs)

    def zlabel(self, title, *args, **kwargs):
        return plt.zlabel(title, *args, **kwargs)

    def xlim(self, *args, **kwargs):
        return plt.xlim(*args, **kwargs)

    def ylim(self, *args, **kwargs):
        return plt.ylim(*args, **kwargs)

    def zlim(self, *args, **kwargs):
        return plt.zlim(*args, **kwargs)

    def legend(self, *args, **kwargs):
        return plt.legend(*args, **kwargs)

    def update_mpl_plot(self):
        plt.gcf().canvas.draw()
        plt.show(block=False)


class MplPlotUtils2D(PlotUtils2D, MplPlotUtilsBase):
    def __init__(self):
        pass

    def init_cr3bp_plot(self, **kwargs):
        plt.figure()
        plt.axis('equal')
        self.default_axes(**kwargs)

    def show(self, **kwargs):
        plt.show(block=False)

    def plot_traj_list_projection(self,
                                  traj_list,
                                  var1='x',
                                  var2='y',
                                  labels=None,
                                  n_highlighted_orbits=None,
                                  cmap_default=None,
                                  cmap_name=None,
                                  cmap_fn=None,
                                  cmap_title=None,
                                  **kwargs):
        lines = []
        step = None

        if cmap_default:
            cmap_name_default = 'jet'

            if cmap_default == 'dy0':
                cmap_fn_default = lambda traj: traj.start().dy
                cmap_title_default = 'Initial dy [nd]'
            elif cmap_default == 'c':
                cmap_fn_default = lambda traj: traj.start().c()
                cmap_title_default = 'Jacobi Constant (C)'
            elif cmap_default == 'period':
                cmap_fn_default = lambda traj: traj.times[1]
                cmap_title_default = 'Period [nd]'
            else:
                raise ValueError('Improper Default Input value %s' % cmap_default)

            cmap_name = cmap_name or cmap_name_default
            cmap_fn = cmap_fn or cmap_fn_default
            cmap_title = cmap_title or cmap_title_default

        scalarMap = None

        if cmap_name and cmap_fn:
            vals = [cmap_fn(traj) for traj in traj_list]
            cmap = plt.get_cmap(cmap_name)
            cNorm = colors.Normalize(vmin=min(vals), vmax=max(vals))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
            scalarMap.set_array(vals)

        if n_highlighted_orbits:
            step = len(traj_list) / n_highlighted_orbits

        for i, traj in enumerate(traj_list):
            if labels:
                kwargs['label'] = labels[i]

            c = None

            if step:
                c = 'b' if i % step == 0 or i == len(traj_list) - 1 else '0.85'

            if c:
                kwargs['c'] = c

            if scalarMap:
                kwargs['c'] = scalarMap.to_rgba(vals[i])

            if i > 0 and 'label' in kwargs:
                kwargs.pop('label')
                traj.name = None

            if 'ls' not in kwargs:
                kwargs['ls'] = 'solid'

            tl = self.plot_traj_projection(traj, var1, var2, show_labels=False, **kwargs)
            lines.append(tl)

        if cmap_title and scalarMap:
            cb = plt.colorbar(scalarMap)
            if cmap_title:
                cb.set_label(cmap_title)

        return lines

    def plot(self, x_vals, y_vals, **kwargs):
        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 0.7
        if 'ls' not in kwargs:
            kwargs['ls'] = 'solid'

        return plt.plot(x_vals, y_vals, **kwargs)[0]

    def plot_li(self, system, li, **kwargs):
        fig = plt.gcf()
        ax = fig.gca()

        if 'marker' not in kwargs:
            kwargs['marker'] = '^'

        if 'c' not in kwargs:
            kwargs['c'] = 'maroon'

        if 'ms' not in kwargs and 'markersize' not in kwargs:
            kwargs['ms'] = 3.0

        kwargs['fillstyle'] = 'none'
        kwargs['lw'] = 0.2
        pos = getattr(system, li.lower()).pos()
        plt.plot(pos[0], pos[1], **kwargs)

    def plot_circle(self, x, y, radius, **kwargs):
        if 'color' not in kwargs and 'c' not in kwargs:
            kwargs['color'] = 'k'

        circle_primary = plt.Circle((x, y), radius, **kwargs)
        plt.gca().add_artist(circle_primary)

    def plot_poincare(self, x, y, **kwargs):
        if 's' not in kwargs:
            kwargs['s'] = 0.15

        return plt.scatter(x, y, **kwargs)

    def plot_jacobi_list(self, system, c_list):
        plt.scatter(range(len(c_list)), c_list)
        plt.xlabel('IC Index')
        plt.ylabel('Jacobi Constant (C)')
        plt.axhline(y=system.l1.c(), ls='dashed', label='C_L1', c='r')
        plt.axhline(y=system.l2.c(), ls='dashed', label='C_L2', c='g')
        plt.axhline(y=system.l3.c(), ls='dashed', label='C_L3', c='c')
        plt.axhline(y=system.l4.c(), ls='dashed', label='C_L4/L5', c='k')
        plt.legend()

    def plot_jacobi_for_states(self, system, states):
        c_list = [s.c() for s in states]
        self.plot_jacobi_list(system, c_list)

    def plot_jacobi_for_arcs(self, system, arcs):
        c_list = [cr3bp.State(system, a.start_state).c() for a in arcs]
        self.plot_jacobi_list(system, c_list)

    def update_plot(self):
        self.update_mpl_plot()

    def plot_accessible_arc_list_w_colormap(self, arc_mngr, q_manager, system, arc, base, time, rank, digraph_args={},
                                            **kwargs):
        accessible_arcs = graph_utils.get_accessible_region_mc(arc_mngr,
                                                               system,
                                                               arc,
                                                               base,
                                                               time,
                                                               **digraph_args)

        terminal_arcs = [a for a in accessible_arcs if a.terminal_arc]
        if terminal_arcs:
            new_arc = terminal_arcs[0]
            self.plot_arc_list(system, [new_arc], c='k')
            return 1000, new_arc

        act_values = q_manager.predict_batch(arc, accessible_arcs, base)
        traj_list_with_q = []
        for a in act_values:
            traj = system.prop(a[1].start_state, 0, a[1].time, n_points=1000)
            traj.metadata['q'] = a[0]
            traj_list_with_q.append(traj)

        cmap_fn = lambda traj: traj.metadata['q']
        sorted(traj_list_with_q, key=cmap_fn, reverse=True)

        if 'cmap_fn' not in kwargs:
            kwargs['cmap_fn'] = cmap_fn

        if 'cmap_name' not in kwargs:
            kwargs['cmap_name'] = 'cool'

            self.plot_traj_list_projection(traj_list_with_q, **kwargs)

        q, new_arc = sorted(act_values, key=lambda av: av[0], reverse=True)[rank]
        return q, new_arc

    def xlabel(self, xlabel, **kwargs):
        self.mpl_xlabel(xlabel, **kwargs)

    def ylabel(self, ylabel, **kwargs):
        self.mpl_ylabel(ylabel, **kwargs)

    def title(self, title, **kwargs):
        self.mpl_title(title, **kwargs)


class MplPlotUtils3D(PlotUtils3D, MplPlotUtilsBase):
    def set_3d_scaling(self):
        ax = plt.gca()

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    def init_cr3bp_plot(self):
        fig = plt.gcf()

        ax = fig.gca(projection='3d')
        ax.set_aspect("equal")

        # ax.set_facecolor('white')

        # I want to update the grid when they add suppose for custom 3d grids check back later on
        # https://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.grid
        ax.grid(True)  # TODO

        # Some configurations found here
        # https://dawes.wordpress.com/2014/06/27/publication-ready-3d-figures-from-matplotlib/

        # White Background
        # ax.xaxis.pane.set_edgecolor('black')
        # ax.yaxis.pane.set_edgecolor('black')
        # ax.zaxis.pane.set_edgecolor('black')

        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        # Font Type
        plt.rc('font', size=9)
        plt.rc('font', family='serif')
        plt.rc('axes', labelsize=10)

        # Tick Label Placement
        [t.set_va('center') for t in ax.get_yticklabels()]
        [t.set_ha('left') for t in ax.get_yticklabels()]
        [t.set_va('center') for t in ax.get_xticklabels()]
        [t.set_ha('right') for t in ax.get_xticklabels()]
        [t.set_va('center') for t in ax.get_zticklabels()]
        [t.set_ha('left') for t in ax.get_zticklabels()]

        # Tick Placement
        ax.xaxis._axinfo['tick']['inward_factor'] = 0
        ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
        ax.yaxis._axinfo['tick']['inward_factor'] = 0
        ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
        ax.zaxis._axinfo['tick']['inward_factor'] = 0
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.4
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.4

        # Not sure if this is necessary
        # ax.xaxis.set_major_locator(MultipleLocator(0.1))
        # ax.yaxis.set_major_locator(MultipleLocator(0.1))
        # ax.zaxis.set_major_locator(MultipleLocator(0.1))

    def plot_sphere(self, x, y, z, radius, **kwargs):
        self.plot_sphere_matplotlib(x, y, z, radius, **kwargs)

    def plot3(self, x_vals, y_vals, z_vals, **kwargs):
        ax = plt.gca()
        return ax.plot(x_vals, y_vals, z_vals, **kwargs)

    def scatter3(self, x_vals, y_vals, z_vals, **kwargs):
        if 's' not in kwargs:
            s = 0.2
        return plt.scatter(x_vals, y_vals, z_vals, **kwargs)

    def default_axes(self, state=None):
        unit = 'nd' if not state or not state.dimensional else 'km'
        fmt_str = '\n\n$%s \; [%s]$'

        ax = plt.gca()
        ax.set_xlabel(fmt_str % ('x', unit), size='large')
        ax.set_ylabel(fmt_str % ('y', unit), size='large')

        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel(fmt_str % ('z', unit), size='large', rotation=90)

    def show(self):
        self.set_3d_scaling()
        plt.show(block=False)

    def update_plot(self):
        self.update_mpl_plot()


# ========================================================================= #
#                             Jacobi Plotting                               #
# ========================================================================= #

def plot_jacobi_error_over_traj(traj, labels=True):
    plt.plot(traj.times, traj.c_error(all_vals=True))
    if labels:
        plt.title("Jacobi Constant Error over Trajectory")
        plt.xlabel('$t \; [nd]$')
        plt.ylabel('$C \;$ Error')


# ========================================================================= #
#                               Tile Plotting                               #
# ========================================================================= #

def tile_graph_grid(edge_length, x_freq=2, y_freq=2, major_freq=2):
    ax = plt.gca()

    x_rounded_lim = [utils.round_to_base(x, edge_length) - edge_length / 2. for x in ax.get_xlim()]
    y_rounded_lim = [utils.round_to_base(y, edge_length) - edge_length / 2. for y in ax.get_ylim()]

    print('X Rounded Lim: {} : {}'.format(x_rounded_lim[0], x_rounded_lim[1]))
    print('Y Rounded Lim: {} : {}'.format(y_rounded_lim[0], y_rounded_lim[1]))

    # x_rounded_lim = [utils.round_to_base(x, edge_length) for x in [-2., 2.]]
    # y_rounded_lim = [utils.round_to_base(y, edge_length) for y in [-2., 2.]]

    x_minor_ticks = np.arange(x_rounded_lim[0], x_rounded_lim[1], step=edge_length)
    x_major_ticks = np.arange(x_rounded_lim[0], x_rounded_lim[1], step=edge_length * major_freq)
    y_minor_ticks = np.arange(y_rounded_lim[0], y_rounded_lim[1], step=edge_length)
    y_major_ticks = np.arange(y_rounded_lim[0], y_rounded_lim[1], step=edge_length * major_freq)

    ax.set_xticks(x_major_ticks)
    ax.set_xticks(x_minor_ticks, minor=True)
    ax.set_yticks(y_major_ticks)
    ax.set_yticks(y_minor_ticks, minor=True)

    ax.grid(which='both', axis='both', linestyle='-')

    # ax.xaxis.set_major_locator(plticker.MultipleLocator(base=edge_length))
    # ax.yaxis.set_major_locator(plticker.MultipleLocator(base=edge_length))
    # ax.grid(which='major', axis='both', linestyle='-')

    for i, x_label in enumerate(ax.get_xticklabels()):
        x_label.set_visible(i % x_freq == 0)
    for j, y_label in enumerate(ax.get_yticklabels()):
        y_label.set_visible(j % y_freq == 0)

    def on_resize(event):
        tile_graph_grid(edge_length, x_freq, y_freq)

    def on_key(event):
        if event.key == 'g':
            tile_graph_grid(edge_length, x_freq, y_freq)

    plt.gcf().canvas.mpl_connect('key_release_event', on_key)


def tile_colormap(edges, name='autumn'):
    colormap_fn = lambda edge: edge.count
    vals = np.array([colormap_fn(edge) for edge in edges])
    c_norm = colors.Normalize(vmin=min(vals), vmax=max(vals))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=plt.get_cmap(name))
    scalar_map.set_array(vals)
    return scalar_map, vals


def tile_patch(var1, var2, base, **kwargs):
    if 'color' not in kwargs:
        kwargs['color'] = 'k'
    if 'c' in kwargs:
        kwargs['color'] = kwargs['c']
        kwargs.pop('c')
    ax = plt.gca()
    ax.add_patch(
        patches.Rectangle(
            (var1 - base / 2, var2 - base / 2),
            base,
            base,
            **kwargs
        )
    )


def tile_patch_from_state_array(st, base, **kwargs):
    tile_patch(utils.round_to_base(st[0], base),
               utils.round_to_base(st[1], base),
               base,
               **kwargs)


def tile_patch_from_state(st, base, **kwargs):
    tile_patch_from_state_array(st.pos(), base, **kwargs)


def tile_patch_from_node(node, base, **kwargs):
    tile_patch(node[0], node[1], base, **kwargs)


def tile_graph_for_pt_from_db(start_point,
                              base,
                              time,
                              x_grid_base=2,
                              y_grid_base=2,
                              colormap_name='autumn',
                              C=None,
                              zvc=False,
                              cache_zvc_path=None):
    system = start_point.system

    start_region_loc = [utils.round_to_base(start_point.x, base), utils.round_to_base(start_point.y, base)]
    edges = arc_db.ArcDb.find_dt_arcs_by_start_state(start_point, base, time)
    traj_list = []

    if not edges:
        print('No edges found for query!')
        return

    for edge in edges:
        traj = system.prop(edge.start_state, 0, float(edge.time))
        traj_list.append(traj)

    # Init Plot
    ax = plt.gca()
    init_cr3bp_plot()

    # Color Map
    scalar_map, color_vals = tile_colormap(edges, name=colormap_name)
    cb = plt.colorbar(scalar_map)
    cb.set_label('Number of Connections')

    for i, traj in enumerate(traj_list):
        c = scalar_map.to_rgba(color_vals[i])
        tile_patch_from_state(traj.end(), base, c=c)

    tile_patch_from_node(start_region_loc, base, c='k')

    if C and zvc:
        args = {
            'default_step': 1e-3
        }

        if cache_zvc_path:
            args['cache_zvc_path'] = cache_zvc_path

        system.zvc(C, **args)

    tile_graph_grid(base, x_grid_base, y_grid_base)
    return traj_list, edges

# ========================================================================= #
#                                 RL Arcs                                   #
# ========================================================================= #
