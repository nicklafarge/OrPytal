def get_plot_utils(engine, planar=False):
    if engine == 'matplotlib':
        from orpytal.plotting import mpl_plotting
        if planar:
            return mpl_plotting.MplPlotUtils2D()
        else:
            return mpl_plotting.MplPlotUtils3D()
    elif engine == 'mayavi':
        from orpytal.plotting.mayavi_plotting import MayaviPlotUtils3D
        return MayaviPlotUtils3D()
    elif engine == 'matlab':
        from orpytal.plotting.matlab_plotting import MatlabPlotUtils3D
        return MatlabPlotUtils3D()
    elif engine == 'plotly':
        from orpytal.plotting import plotly_plotting
        if planar:
            return plotly_plotting.PlotlyPlotUtils2D()
        else:
            return plotly_plotting.PlotlyPlotUtils3D()
    else:
        raise NameError('Invalid Plot Library. Valid keys: matplotlib, mayavi, matlab, plotly')