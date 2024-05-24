import warnings
import numpy as np
from ctypes import c_int, c_float, c_void_p, POINTER

from audioflux.fftlib import get_fft_lib

__all__ = [
    'ascontiguous_T',
    'ascontiguous_swapaxex',
    'format_channel',
    'revoke_channel',
    'check_audio',
    'check_audio_length',
    'synth_f0'
]


def ascontiguous_T(X, dtype=None, *args, **kwargs):
    """
    The transposed array, and then return a contiguous array in memory.

    Parameters
    ----------
    X: ndarray
        Input array.
    dtype: str or dtype object, optional
        Data-type of returned array.

    Returns
    -------
    out: ndarray
        return a transposed and contiguous array
    """
    return np.ascontiguousarray(X.T, dtype=dtype, *args, **kwargs)


def ascontiguous_swapaxex(X, axis1, axis2, dtype=None, *args, **kwargs):
    """
    Interchange two axes of an array, and then return a contiguous array in memory.

    Parameters
    ----------
    X: ndarray
        Input array.

    dtype: str or dtype object, optional
        Data-type of returned array.

    Returns
    -------
    out: ndarray
        return a transposed and contiguous array
    """
    return np.ascontiguousarray(np.swapaxes(X, axis1, axis2), dtype=dtype, *args, **kwargs)


def format_channel(X, last_fixed_ndim):
    shape = X.shape
    return X.reshape((-1, *shape[-last_fixed_ndim:])), shape[:-last_fixed_ndim]


def revoke_channel(X, target_channel_shape, last_fixed_ndim):
    return X.reshape((*target_channel_shape, *X.shape[-last_fixed_ndim:]))


def check_audio(X, is_mono=True):
    """
    check audio data

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio array

    is_mono: bool
        is mono

    Returns
    -------
    out: bool
    """
    if not isinstance(X, np.ndarray):
        raise ValueError(f'Audio data must be of type np.ndarray')

    if not X.flags.c_contiguous:
        raise ValueError("Audio data must be C-contiguous")

    if X.ndim == 0:
        raise ValueError(f'Audio data must have at least one dimension')

    if is_mono and X.ndim != 1:
        raise ValueError(f'Audio data must be a monophonic, there is only 1 dimension'
                         f', given X.ndim={X.ndim}, X.shape={X.shape}')

    return True


def check_audio_length(X, radix2_exp):
    data_len = X.shape[-1]
    fft_length = 1 << radix2_exp
    if data_len < fft_length:
        pad_len = fft_length - data_len
        warnings.warn(f'The audio length={data_len} is not enough for fft_length={fft_length}(2**radix2_exp), '
                      f'and {pad_len} zeros are automatically filled after the audio')
        X = np.pad(X, (*[(0, 0)] * (X.ndim - 1), (0, pad_len)))
    elif data_len > fft_length:
        warnings.warn(f'fft_length={fft_length}(2**radix2_exp) is too small for data_arr length={data_len}, '
                      f'only the first fft_length={fft_length} data are valid')
        X = X[..., :fft_length].copy()
    return X


def synth_f0(times, frequencies, samplate, amplitudes=None):
    """
    Generate an audio array based on the frequency f0.

    Parameters
    ----------
    times: ndarray  [shape=(n)]
        Time points for each frequency, in seconds.

    frequencies: ndarray [shape=(n)]
        Array of frequencies, in Hz.

    samplate: int
        The output sampling rate.

    amplitudes: ndarray [shape=(n)]
        The amplitude of each frequency, ranging from 0 to 1.

        Default is None, which means that the amplitude is 1. Like: `np.ones((n,))`

    Returns
    -------
    out: ndarray
        Return the audio array generated based on the frequencies


    Examples
    --------

    >>> import audioflux as af
    >>> import numpy as np
    >>> f0_arr = np.ones((1024,)) * 220
    >>> times = np.arange(0, f0_arr.shape[0]) * (1024 / 32000)
    >>> amplitude_arr = np.ones_like(f0_arr) * 0.4
    >>> audio_arr = af.utils.synth_f0(times, f0_arr, 32000, amplitude_arr)

    """
    times = np.asarray(times, dtype=np.float32, order='C')
    if times.ndim != 1:
        raise ValueError(f"times[ndim={times.ndim}] must be a 1D array")

    frequencies = np.asarray(frequencies, dtype=np.float32, order='C')
    if frequencies.ndim != 1:
        raise ValueError(f"frequencies[ndim={frequencies.ndim}] must be a 1D array")

    if times.shape[0] != frequencies.shape[0]:
        raise ValueError(f"The lengths of times and frequencies must be the same.")

    if amplitudes is not None:
        amplitudes = np.asarray(amplitudes, dtype=np.float32, order='C')
        if amplitudes.ndim != 1:
            raise ValueError(f"amplitudes[ndim={amplitudes.ndim}] must be a 1D array")

        if amplitudes.shape[0] != frequencies.shape[0]:
            raise ValueError(f"The lengths of amplitudes and frequencies must be the same.")

    fn = get_fft_lib()['util_synthF0']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int, c_int,
        POINTER(c_void_p) if amplitudes is None else np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    ]
    fn.restype = c_void_p

    length = times.shape[0]
    length1 = round(np.max(times) * samplate)

    p = fn(times, frequencies, c_int(length), c_int(samplate), amplitudes)
    ret = np.frombuffer((c_float * length1).from_address(p), np.float32).copy()
    return ret