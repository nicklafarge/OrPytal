def get_plot_utils(engine, planar=False):
    if engine == 'matplotlib':
        import mpl_plotting
        if planar:
            return mpl_plotting.MplPlotUtils2D()
        else:
            return mpl_plotting.MplPlotUtils3D()
    elif engine == 'matlab':
        from matlab_plotting import MatlabPlotUtils3D
        return MatlabPlotUtils3D()
    elif engine == 'plotly':
        import plotly_plotting
        if planar:
            return plotly_plotting.PlotlyPlotUtils2D()
        else:
            return plotly_plotting.PlotlyPlotUtils3D()
    else:
        raise NameError('Invalid Plot Library. Valid keys: matplotlib, matlab, plotly')