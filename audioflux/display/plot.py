import warnings
import matplotlib.pyplot as plt
import numpy as np

from .display import fill_spec, fill_plot, fill_wave

__all__ = ['Plot']


class Plot(object):
    """
    Create a set of subplots for audioflux.

    Parameters
    ----------
    nrows: int, default: 1
        Number of rows of the subplot grid.

    ncols: int, default: 1
        Number of columns of the subplot grid.

    sharex  bool or {'none', 'all', 'row', 'col'}, default: False
        Controls sharing of properties among x (*sharex*)

    sharey: bool or {'none', 'all', 'row', 'col'}, default: False
        Controls sharing of properties among y (*sharey*)

    fig_width: float
        Width in inches.

    fig_height: float
        Height in inches.

    fig_kw:
        Finger params

    """

    def __init__(self, nrows=1, ncols=1, sharex=False, sharey=False, fig_width=8, fig_height=8, fig_kw=None):
        self.nrows = nrows
        self.ncols = ncols
        fig_kw = fig_kw or {}
        fig_kw['figsize'] = (fig_width, fig_height)
        self.fig, self.ax = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey, **fig_kw)

    def get_axes(self, row, col):
        """
        Get axes

        Parameters
        ----------
        row: int
            Index of row of the subplot grid.

        col: int
            Index of col of the subplot grid.

        Returns
        -------
        ax: matplotlib.axes.Axes
        """
        if self.nrows == 1 and self.ncols == 1:
            return self.ax
        elif self.nrows == 1:
            return self.ax[col]
        elif self.ncols == 1:
            return self.ax[row]
        return self.ax[row][col]

    def add_spec_data(self, data, x_coords=None, y_coords=None, scale='linear', row_idx=0, col_idx=0, title='',
                      *, show_colorbar=True, axis_option='on'):
        """
        Add spectrogram data

        Parameters
        ----------
        data: np.ndarray [shape=(fre, time)]
            The matrix of spectrogram

        x_coords, y_coords: None or np.ndarray [shape=(n,)]
            Array of X-axis/Y-axis ticks

            If not provided, it will be divided equally according to the size of the matrix

        scale: str
            The scale of Y-axis

            - 'linear': the spectrum is displayed on a linear scale.
            - 'log': the spectrum is displayed on a log scale.

            If ``y_coords`` is None, it will always be linear

        row_idx:
            Index of row of the subplot grid.

        col_idx:
            Index of col of the subplot grid.

        title: str
            Set a title for the Axes.

        show_colorbar: bool
            Whether to display the colorbar

        axis_option: bool or str
            See `matplotlib._AxesBase.axis()` for details.

        Returns
        -------
        ax: matplotlib.axes.Axes
        """
        if x_coords is None:
            x_axis = None
        else:
            x_axis = 'time'
        if y_coords is None and scale != 'linear':
            scale = 'linear'
            warnings.warn(f'If `y_coords` is None, `scale` must be linear')
        ax = self.get_axes(row_idx, col_idx)
        img = fill_spec(data, ax, x_coords=x_coords, y_coords=y_coords, x_axis=x_axis, y_axis=scale, title=title)
        if show_colorbar:
            self.fig.colorbar(img, ax=ax)
        plt.axis(axis_option)
        return ax

    def add_wave_data(self, data, samplate=32000, row_idx=0, col_idx=0):
        """
        Add wave data

        Parameters
        ----------
        data: np.ndarray [shape=(n,)]
            Input audio data

        samplate: int
            Sampling rate of the incoming audio

        row_idx:
            Index of row of the subplot grid

        col_idx:
            Index of col of the subplot grid


        Returns
        -------
        ax: matplotlib.axes.Axes
        """
        axes = self.get_axes(row=row_idx, col=col_idx)
        return fill_wave(data, samplate=samplate, axes=axes)

    def add_plot(self, x, y, label='', row_idx=0, col_idx=0, is_legend=True, *,
                 x_lims=None, y_lims=None, y_blank_threshold=0.15):
        """
        Add plot data

        Parameters
        ----------
        x: np.ndarray [shape=(time,)]
            The horizontal(time) coordinates of the data points.

        y: np.ndarray [shape=(data,)]
            The vertical coordinates of the data points.

        label: str
            Set a label for the Axes.

        row_idx, col_idx: int
            Index of row/col of the subplot grid.

        x_lims: None or (left: float, right: float)
            Set the x-axis view limits.

            * None: Automatically get the minimum and maximum values for the left and right of the x-axis coordinates.
            * (left: float, right: float): The left and right of the x-axis coordinates.

        y_lims: None or (bottom: float, top: float)
            Set the y-axis view limits.

            * None: Automatically get the minimum and maximum values for the bottom and top of the x-axis coordinates.
            * (bottom: float, top: float): The bottom and top of the y-axis coordinates.

        y_blank_threshold: float
            The space between the minimum and maximum values of the data and the y-axis scale.

        Returns
        -------
        ax: matplotlib.axes.Axes
        """
        axes = self.get_axes(row_idx, col_idx)

        axes = fill_plot(x, y, axes=axes, label=label, is_legend=is_legend,
                         x_lims=x_lims, y_lims=y_lims, y_blank_threshold=y_blank_threshold)
        return axes

    def show(self):
        """
        Show plot
        """
        plt.show()

    def save(self, fname, *, dpi='figure', format=None, metadata=None,
             bbox_inches=None, pad_inches=0.1,
             facecolor='auto', edgecolor='auto',
             backend=None, **kwargs
             ):
        """
        Save plot to file

        See `matplotlib.savefig`
        """

        plt.savefig(fname, dpi=dpi, format=format, metadata=metadata,
                    bbox_inches=bbox_inches, pad_inches=pad_inches,
                    facecolor=facecolor, edgecolor=edgecolor,
                    backend=backend, **kwargs
                    )

    def close(self, fig='all'):
        """
        Close a figure window.

        Parameters
        ----------
        fig : None or int or str or `.Figure`
            The figure to close. There are a number of ways to specify this:

            - *None*: the current figure
            - `.Figure`: the given `.Figure` instance
            - ``int``: a figure number
            - ``str``: a figure name
            - 'all': all figures
        """
        plt.close(fig)
