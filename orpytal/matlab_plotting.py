########### Standard ###########
import logging

########### Local ###########
import orpytal.orbit
import orpytal.planet_constants
from orpytal.common import units

########### External ###########
import matlab.engine

logging.info("Starting Matlab Engine....")
eng = matlab.engine.start_matlab()
logging.info("Matlab Engine Started!")


def eval(cmd, nargout=0):
    eng.eval(cmd, nargout=nargout)


def figure(fig_number=None):
    if fig_number:
        eng.figure(float(fig_number))
    else:
        eng.figure()


def init_cr3bp_3d_plot(fig_number=None):
    figure(fig_number)

    eval("hold on;")
    eval("axis equal;")
    default_3d_axes()


def plot_sphere(x, y, z, radius):
    if isinstance(radius, units.Quantity):
        radius = radius.to('km').m
    eng.workspace['radius'] = radius
    eng.workspace['x'] = float(x)
    eng.workspace['y'] = float(y)
    eng.workspace['z'] = float(z)
    eval("[sx sy sz] = sphere;")
    eval("surf(sx*radius + x, sy*radius + y, sz*radius + z, 'FaceColor', [0 0 0]);")


def plot_primary(primary):
    if isinstance(primary, orbit.Orbit):
        radius = primary.central_body.radius
    elif isinstance(primary, planet_constants.CentralBody):
        radius = primary.radius
    else:
        radius = primary
    return plot_sphere(0, 0, 0, radius)


def plot_p1(system):
    return plot_sphere(-system.mu, 0, 0, system.primary_body.radius / system.lstar)


def legend(*args):
    if not args:
        eng.legend()
    else:
        eng.legend(*args)


def _prepare_plot(x_vals, y_vals, z_vals, **kwargs):
    x_vals = _to_list(x_vals)
    y_vals = _to_list(y_vals)
    z_vals = _to_list(z_vals)

    x_vals = matlab.double(x_vals)
    y_vals = matlab.double(y_vals)
    z_vals = matlab.double(z_vals)
    args = format_args(**kwargs)
    return x_vals, y_vals, z_vals, args


def _to_list(val):
    if not hasattr(val, '__len__'):
        val = [val]
    return val


def scatter3(x_vals, y_vals, z_vals, **kwargs):
    x_vals, y_vals, z_vals, args = _prepare_plot(x_vals, y_vals, z_vals, **kwargs)

    args = replace_arg_key(args, 'Color', 'MarkerEdgeColor')

    return eng.scatter3(x_vals, y_vals, z_vals, *args, nargout=1)


def replace_arg_key(args, key, new_key):
    args_list = list(args)
    for i, arg in enumerate(args_list):
        if arg == key:
            args_list[i] = new_key
    return tuple(args_list)


def format_args(**kwargs):
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


def plot3(x_vals, y_vals, z_vals, **kwargs):
    x_vals, y_vals, z_vals, args = _prepare_plot(x_vals, y_vals, z_vals, **kwargs)
    return eng.plot3(x_vals, y_vals, z_vals, *args, nargout=1)


def scatter3_states(states, **kwargs):
    return scatter3(
        [st.x for st in states],
        [st.y for st in states],
        [st.z for st in states],
        **kwargs
    )


def plot_3d_traj_list(traj_list, **kwargs):
    logging.debug("Plotting list of 3D Trajectories....")
    for i, traj in enumerate(traj_list):
        if (i + 1) % 10 == 0:
            logging.info("Traj {} / {}".format(i + 1, len(traj_list)))
        plot_3d_traj(traj, **kwargs)
    logging.debug("Plotting list of 3D Trajectories done!")


def plot_3d_traj(traj, **kwargs):
    if 'color' in traj.metadata and 'c' not in kwargs:
        kwargs['c'] = traj.metadata['color']

    xs, ys, zs = traj.inertial()
    return plot3(xs, ys, zs, **kwargs)


def title(title, latex=True):
    args = ()
    if latex:
        args = args + ("Interpreter", 'Latex')
    return eng.title(title, *args)


def xlabel(xlabel, latex=True):
    args = ()
    if latex:
        args = args + ("Interpreter", 'Latex')

    return eng.xlabel(xlabel, *args)


def ylabel(ylabel, latex=True):
    args = ()
    if latex:
        args = args + ("Interpreter", 'Latex')

    return eng.ylabel(ylabel, *args)


def zlabel(zlabel, latex=True):
    args = ()
    if latex:
        args = args + ("Interpreter", 'Latex')

    return eng.zlabel(zlabel, *args)


def default_3d_axes():
    unit = 'km'
    fmt_str = '$\hat{%s}$ [%s]'

    xlabel(fmt_str % ('x', unit), latex=True)
    ylabel(fmt_str % ('y', unit), latex=True)
    zlabel(fmt_str % ('z', unit), latex=True)


def grid():
    eng.eval("grid on;", nargout=0)


def show():
    pass
