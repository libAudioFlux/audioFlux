from ctypes import c_int, c_float, POINTER, pointer
import numpy as np

from audioflux.fftlib import get_fft_lib

__all__ = [
    'power_to_db',
    'power_to_abs_db',
    'mag_to_abs_db',
    'log_compress',
    'log10_compress',
    'delta',
    'get_phase',
]


def power_to_db(X, min_db=-80):
    """
    Convert power spectrogram to relative decibel(dB) spectrogram.

    Parameters
    ----------
    X: np.ndarray
        Input array

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray
        relative dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_powerToDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    length = arr.size
    db_arr = np.zeros(arr.shape, dtype=np.float32)
    fn(arr, c_int(length), c_float(min_db), db_arr)
    return db_arr.reshape(shape)


def power_to_abs_db(X, fft_length=4096, is_norm=False, min_db=-80):
    """
    Convert power spectrogram to absolute decibel(dB) spectrogram

    Parameters
    ----------
    X: np.ndarray
        Input array

    fft_length: int
        fft length

    is_norm: bool
        Whether to use normalization

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray
        absolute dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_powerToAbsDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    length = arr.size
    db_arr = np.zeros(arr.shape, dtype=np.float32)
    fn(arr, c_int(length), c_int(fft_length), c_int(int(is_norm)), c_float(min_db), db_arr)
    return db_arr.reshape(shape)


def mag_to_abs_db(X, fft_length=4096, is_norm=False, min_db=-80):
    """
    Convert magnitude spectrogram to absolute decibel(dB) spectrogram

    Parameters
    ----------
    X: np.ndarray
        Input array

    fft_length: int
        fft length

    is_norm: bool
        Whether to use normalization

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray
        absolute dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_magToAbsDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    length = arr.size
    db_arr = np.zeros(arr.shape, dtype=np.float32)
    fn(arr, c_int(length), c_int(fft_length), c_int(int(is_norm)), c_float(min_db), db_arr)
    return db_arr.reshape(shape)


def log_compress(X, gamma=1.0):
    """
    log compression

    Parameters
    ----------
    X: np.ndarray [shape=(fre,) or (fre, time)]
        Input array

    gamma: float
        1/10/20/...

    Returns
    -------
    out: np.ndarray [shape=(fre,) or (fre, time)]
        log compressed array
    """
    X = np.asarray(X, dtype=np.float32, order='C')

    if X.ndim > 2:
        raise ValueError(f'log10_compress is only defined for 1D arrays')

    fn = get_fft_lib()['util_logCompress']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        POINTER(c_float),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    ret_arr = np.zeros_like(X, dtype=np.float32)
    fn(X, pointer(c_float(gamma)), c_int(X.shape[0]), ret_arr)
    return ret_arr


def log10_compress(X, gamma=1.0):
    """
    log10 compression

    Parameters
    ----------
    X: np.ndarray [shape=(fre,) or (fre, time)]
        Input array

    gamma: float
        1/10/20/...

    Returns
    -------
    out: np.ndarray [shape=(fre,) or (fre, time)]
        log10 compressed array
    """
    X = np.asarray(X, dtype=np.float32, order='C')

    if X.ndim > 2:
        raise ValueError(f'log10_compress is only defined for 1D arrays')

    fn = get_fft_lib()['util_log10Compress']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        POINTER(c_float),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    ret_arr = np.zeros_like(X, dtype=np.float32)
    fn(X, pointer(c_float(gamma)), c_int(X.shape[0]), ret_arr)
    return ret_arr


def delta(X, order=9):
    """
    Compute delta features

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        Input 1D array

    order: int
        must odd

    Returns
    -------
    out: np.ndarray
        delta array
    """
    if X.ndim > 1:
        raise ValueError(f'delta is only defined for 1D arrays')

    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_delta']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    length = X.shape[0]
    ret_arr = np.zeros_like(X, dtype=np.float32)
    fn(X, c_int(length), c_int(order), ret_arr)
    return ret_arr


def get_phase(D):
    """
    Extract phase data from a complex-valued spectrogram D.

    Parameters
    ----------
    D: np.ndarray [shape=(fre, time)]
        a complex-valued spectrogram

    Returns
    -------
    out: np.ndarray [shape=(fre, time)]
        phase data
    """
    if not np.iscomplexobj(D):
        raise ValueError(f'D must be type of complex')

    m_real_arr = D.real.copy()
    m_imag_arr = D.imag.copy()
    m_real_arr[m_real_arr < 1e-16] = 1e-16
    m_phase_arr = np.arctan2(m_imag_arr, m_real_arr)
    return m_phase_arr
