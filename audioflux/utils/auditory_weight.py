from ctypes import c_int

import numpy as np
from audioflux.fftlib import get_fft_lib

__all__ = [
    "auditory_weight_a",
    "auditory_weight_b",
    "auditory_weight_c",
    "auditory_weight_d"
]


def auditory_weight_a(fre_arr):
    """
    Compute the weight-a of a set of frequencies.

    Parameters
    ----------
    fre_arr: np.ndarray [shape=(n,)]
        Input frequencies array

    Returns
    -------
    out: np.ndarray [shape=(n,)]
    """
    fre_arr = np.asarray(fre_arr, dtype=np.float32, order='C')
    if fre_arr.ndim != 1:
        raise ValueError(f'auditory_weight_d is only defined for 1D arrays')

    fn = get_fft_lib()['auditory_weightA']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]
    length = fre_arr.shape[0]
    ret = np.zeros_like(fre_arr, dtype=np.float32)
    fn(fre_arr, c_int(length), ret)
    return ret


def auditory_weight_b(fre_arr):
    """
    Compute the weight-b of a set of frequencies.

    Parameters
    ----------
    fre_arr: np.ndarray [shape=(n,)]
        Input frequencies array

    Returns
    -------
    out: np.ndarray [shape=(n,)]
    """
    fre_arr = np.asarray(fre_arr, dtype=np.float32, order='C')
    if fre_arr.ndim != 1:
        raise ValueError(f'auditory_weight_d is only defined for 1D arrays')

    fn = get_fft_lib()['auditory_weightB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]
    length = fre_arr.shape[0]
    ret = np.zeros_like(fre_arr, dtype=np.float32)
    fn(fre_arr, c_int(length), ret)
    return ret


def auditory_weight_c(fre_arr):
    """
    Compute the weight-c of a set of frequencies.

    Parameters
    ----------
    fre_arr: np.ndarray [shape=(n,)]
        Input frequencies array

    Returns
    -------
    out: np.ndarray [shape=(n,)]
    """
    fre_arr = np.asarray(fre_arr, dtype=np.float32, order='C')
    if fre_arr.ndim != 1:
        raise ValueError(f'auditory_weight_d is only defined for 1D arrays')

    fn = get_fft_lib()['auditory_weightC']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]
    length = fre_arr.shape[0]
    ret = np.zeros_like(fre_arr, dtype=np.float32)
    fn(fre_arr, c_int(length), ret)
    return ret


def auditory_weight_d(fre_arr):
    """
    Compute the weight-d of a set of frequencies.

    Parameters
    ----------
    fre_arr: np.ndarray [shape=(n,)]
        Input frequencies array

    Returns
    -------
    out: np.ndarray [shape=(n,)]
    """
    fre_arr = np.asarray(fre_arr, dtype=np.float32, order='C')
    if fre_arr.ndim != 1:
        raise ValueError(f'auditory_weight_d is only defined for 1D arrays')

    fn = get_fft_lib()['auditory_weightD']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]
    length = fre_arr.shape[0]
    ret = np.zeros_like(fre_arr, dtype=np.float32)
    fn(fre_arr, c_int(length), ret)
    return ret
