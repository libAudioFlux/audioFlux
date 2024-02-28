from ctypes import c_float, c_int, POINTER, pointer

from audioflux.fftlib import get_fft_lib

__all__ = [
    'queue_fre2',
    'queue_fre3',
]


def queue_fre2(fre1, fre2):
    """
    Calculate the queue frequency of two frequencies.

    Parameters
    ----------
    fre1: float
        frequency1

    fre2: float
        frequency2

    Returns
    -------
    k1: int
        The ratio1 relationship between fre1 and fre2.
    k2: int
        The ratio2 relationship between fre1 and fre2.
    n: float
        Base frequency of fre1 and fre2.
    """
    fn = get_fft_lib()['__queue_fre2']
    fn.argtypes = [c_float, c_float, POINTER(c_int), POINTER(c_int)]
    fn.restype = c_float

    p_k1_arr = pointer(c_int())
    p_k2_arr = pointer(c_int())
    n = fn(c_float(fre1), c_float(fre2), p_k1_arr, p_k2_arr)
    return p_k1_arr[0], p_k2_arr[0], n


def queue_fre3(fre1, fre2, fre3):
    """
    Calculate the queue frequency of three frequencies.

    Parameters
    ----------
    fre1: float
        frequency1

    fre2: float
        frequency2

    fre3: float
        frequency3

    Returns
    -------
    s1: int
    s2: int
    k1: int
        The ratio1 relationship between fre1, fre2, and fre3.
    k2: int
        The ratio2 relationship between fre1, fre2, and fre3.
    k3: int
        The ratio3 relationship between fre1, fre2, and fre3.
    n: float
        Base frequency of fre1, fre2, and fre3.

    """
    fn = get_fft_lib()['__queue_fre3']
    fn.argtypes = [
        c_float, c_float, c_float,
        POINTER(c_int), POINTER(c_int),
        POINTER(c_int), POINTER(c_int), POINTER(c_int)
    ]
    fn.restype = c_float

    p_s1_arr = pointer(c_int())
    p_s2_arr = pointer(c_int())
    p_k1_arr = pointer(c_int())
    p_k2_arr = pointer(c_int())
    p_k3_arr = pointer(c_int())
    n = fn(c_float(fre1), c_float(fre2), c_float(fre3),
           p_s1_arr, p_s2_arr, p_k1_arr, p_k2_arr, p_k3_arr)
    return p_s1_arr[0], p_s2_arr[0], p_k1_arr[0], p_k2_arr[0], p_k3_arr[0], n
