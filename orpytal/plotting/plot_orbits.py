from orpytal import frames

plot_utils = None


def plot_orbit(orbit, frame=frames.InertialFrame, planar=None):
    # traj = orbit.analytic_propagate_full_orbit()
    traj = orbit.propagate_orbit()
    plot_utils = get_plot_utils(frame=frame)
    plot_utils.init_plot(frame=frame)
    plot_utils.plot_primary(traj)
    plot_utils.plot_traj(traj)
    plot_utils.show()


def get_plot_utils(frame, planar=False):
    from orpytal.plotting import plotly_plotting
    if frame == frames.InertialFrame or not planar:
        return plotly_plotting.PlotlyPlotUtils3D()
    elif frame == frames.PerifocalFrame:
        return plotly_plotting.PlotlyPlotUtils2D()
    else:
        raise ValueError(f"Cannot plot in the {frame.name} with planar={planar}")
