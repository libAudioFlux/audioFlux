import warnings
import numpy as np

__all__ = [
    'ascontiguous_T',
    'ascontiguous_swapaxex',
    'format_channel',
    'revoke_channel',
    'check_audio',
    'check_audio_length'
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
