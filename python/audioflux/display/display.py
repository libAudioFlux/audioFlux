import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as plaxes
from matplotlib.cm import get_cmap
from matplotlib.ticker import Formatter, ScalarFormatter, MaxNLocator, SymmetricalLogLocator, FixedLocator
from audioflux.utils import midi_to_note


def _axis_scale(axes, ax_name, ax_type):
    """Set the x/y-axis scale."""

    if ax_name == "x":
        scaler = axes.set_xscale
    else:
        scaler = axes.set_yscale

    kwargs = {}
    if ax_type == "log":
        mode = "symlog"
        kwargs['base'] = 2
        kwargs['linthresh'] = 64
        kwargs['linscale'] = 0.5
    else:
        mode = "linear"
    scaler(mode, **kwargs)


def _axes_check(axes):
    """Check if "axes" is an instance of an axis object. If not, use `subplot`."""
    if axes is None:
        axes = plt.subplot()
    elif not isinstance(axes, plaxes.Axes):
        raise ValueError("`axes` must be an instance of `matplotlib.axes.Axes`")
    return axes


def _axis_decorate(axis, ax_type, coords):
    """Configure axis tickers, formatters, locators, and labels"""

    if ax_type is None:
        if len(coords) <= 2:
            axis.set_ticks(coords)
        axis.set_label_text("")

    elif ax_type == "time":
        axis.set_major_formatter(TimeFormatter(unit=None, lag=False))
        axis.set_major_locator(MaxNLocator(prune=None, steps=[1, 1.5, 5, 6, 10]))
        axis.set_label_text("Time")

    elif ax_type == "log":
        axis.set_major_formatter(ScalarFormatter())
        axis.set_major_locator(SymmetricalLogLocator(axis.get_transform()))

    elif ax_type == "linear":
        axis.set_major_formatter(ScalarFormatter())

    elif ax_type == "chroma":
        coords_len = len(coords)
        if (coords_len - 1) % 12 != 0:
            raise ValueError(f'The number={coords_len - 1} of y-axis scales of chroma must be a multiple of 12')
        bin_per_tone = (coords_len - 1) // 12
        axis.set_major_formatter(ChromaFormatter(bin_per_tone=bin_per_tone))

        degrees = np.array([0, 2, 4, 5, 7, 9, 11])
        axis.set_major_locator(FixedLocator(degrees * bin_per_tone))
        axis.set_label_text("Pitch class")
    else:
        raise ValueError("Unsupported axis type: {}".format(ax_type))


def _same_axes(x_axis, y_axis, xlim, ylim):
    """Check if two axes are the same, used to determine squared plots"""
    axes_same_and_not_none = (x_axis == y_axis) and (x_axis is not None)
    axes_same_lim = xlim == ylim
    return axes_same_and_not_none and axes_same_lim


class TimeFormatter(Formatter):

    def __init__(self, lag=False, unit=None):

        if unit not in ["s", "ms", None]:
            raise ValueError("Unknown time unit: {}".format(unit))

        self.unit = unit
        self.lag = lag

    def __call__(self, x, pos=None):
        """Return the time format as pos"""

        _, dmax = self.axis.get_data_interval()
        vmin, vmax = self.axis.get_view_interval()

        # In lag-time axes, anything greater than dmax / 2 is negative time
        if self.lag and x >= dmax * 0.5:
            # In lag mode, don't tick past the limits of the data
            if x > dmax:
                return ""
            value = np.abs(x - dmax)
            # Do we need to tweak vmin/vmax here?
            sign = "-"
        else:
            value = x
            sign = ""

        if self.unit == "s":
            s = "{:.3g}".format(value)
        elif self.unit == "ms":
            s = "{:.3g}".format(value * 1000)
        else:
            if vmax - vmin > 3600:
                # Hours viz
                s = "{:d}:{:02d}:{:02d}".format(
                    int(value / 3600.0),
                    int(np.mod(value / 60.0, 60)),
                    int(np.mod(value, 60)),
                )
            elif vmax - vmin > 60:
                # Minutes viz
                s = "{:d}:{:02d}".format(int(value / 60.0), int(np.mod(value, 60)))
            elif vmax - vmin >= 1:
                # Seconds viz
                s = "{:.2g}".format(value)
            else:
                # Milliseconds viz
                s = "{:.3f}".format(value)

        return "{:s}{:s}".format(sign, s)


class ChromaFormatter(Formatter):

    def __init__(self, bin_per_tone=1):
        self.bin_per_tone = bin_per_tone

    def __call__(self, x, pos=None):
        x = int(x) // self.bin_per_tone
        r = midi_to_note(x, is_octave=False)
        return r


def fill_spec(
        data,
        axes=None,
        x_coords=None,
        y_coords=None,
        x_axis=None,
        y_axis=None,
        title=''
):
    """
    Display a spectrogram data

    Parameters
    ----------
    data: np.ndarray [shape=(fre, time)]
        The matrix of spectrogram

    axes: matplotlib.axes.Axes or None
        Axes to plot on instead of the default `plt.subplot()`.

    x_coords, y_coords: None or np.ndarray [shape=(n,)]
        Array of X-axis/Y-axis ticks

        If not provided, it will be divided equally according to the size of the matrix

    x_axis, y_axis: str or None
        The scale of X-axis/Y-axis

        - None: Displayed on a linear scale.
        - 'linear': Displayed on a linear scale.
        - 'log': Displayed on a log scale.
        - 'chroma': Pitches are determined by the chroma filters.
        - 'time': Markers are shown as seconds.

    title: str
        The title of subplot

    Returns
    -------
    matplotlib.collections.QuadMesh
    """
    if np.iscomplexobj(data):
        warnings.warn('Display after performing abs on complex numbers')
        data = np.abs(data)

    if data.ndim != 2:
        raise ValueError(f"data[ndim={data.ndim}] must be a 2D array")

    if x_coords is None:
        x_coords = np.arange(data.shape[-1] + 1)
    if y_coords is None:
        y_coords = np.arange(data.shape[-2] + 1)

    if y_axis == 'chroma':
        y_coords = np.arange(data.shape[-2] + 1)

    cmap = get_cmap('plasma')
    collection = axes.pcolormesh(x_coords, y_coords, data, cmap=cmap)

    axes.set_xlim(np.min(x_coords), np.max(x_coords))
    axes.set_ylim(np.min(y_coords), np.max(y_coords))

    _axis_scale(axes, "x", x_axis)
    _axis_scale(axes, "y", y_axis)

    _axis_decorate(axes.xaxis, x_axis, x_coords)
    _axis_decorate(axes.yaxis, y_axis, y_coords)

    if title:
        axes.set_title(title)

    return collection


def fill_plot(x, y, axes=None, label='', is_legend=True, *,
              x_lims=None, y_lims=None, y_blank_threshold=0.15):
    """
    Display a plot data

    Parameters
    ----------
    x: np.ndarray [shape=(time,)]
        The horizontal(time) coordinates of the data points.

    y: np.ndarray [shape=(data,)]
        The vertical coordinates of the data points.

    axes: matplotlib.axes.Axes or None
        Axes to plot on instead of the default `plt.subplot()`.

    label: str
        Set a label for the Axes.

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
    axes: matplotlib.axes.Axes
    """
    axes = _axes_check(axes)
    if x.ndim != 1:
        raise ValueError(f"x[ndim={x.ndim}] must be a 1D array")
    if y.ndim != 1:
        raise ValueError(f"y[ndim={y.ndim}] must be a 1D array")

    if not x_lims:
        x_lims = (np.min(x), np.max(x))
    x_lims = tuple(x_lims)

    if not y_lims:
        y_min = np.min(y)
        y_max = np.max(y)
        y_blank_len = np.abs(y_max - y_min) * y_blank_threshold
        y_lims = (y_min - y_blank_len, y_max + y_blank_len)
    y_lims = tuple(y_lims)

    axes.set_xlim(*x_lims)
    axes.set_ylim(*y_lims)

    axes.plot(x, y, label=label)
    if is_legend:
        axes.legend()
    return axes


def fill_wave(data, samplate=32000, axes=None):
    """
    Display a wave data

    Parameters
    ----------
    data: np.ndarray [shape=(n,)]
        Input audio data

    axes: matplotlib.axes.Axes or None
        Axes to plot on instead of the default `plt.subplot()`.

    samplate: int
        Sampling rate of the incoming audio

    Returns
    -------
    axes: matplotlib.axes.Axes
    """
    if data.ndim != 1:
        raise ValueError(f"data[ndim={data.ndim}] must be a 1D array")

    times = np.arange(data.shape[-1]) / samplate
    return fill_plot(times, data, axes=axes,
                     x_lims=(times[0], times[-1]),
                     is_legend=False, y_blank_threshold=0.15)
