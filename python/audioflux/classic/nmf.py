from ctypes import c_int, POINTER, c_float, pointer

import numpy as np

from audioflux.fftlib import get_fft_lib

__all__ = [
    'nmf',
]


def nmf(X, k, max_iter=300, tp=0, thresh=1e-3, norm=0):
    X = np.asarray(X, dtype=np.float32, order='C')
    if X.ndim != 2:
        raise ValueError(f"X[ndim={X.ndim}] must be a 2D array")

    fn = get_fft_lib()['nmf']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
        c_int, c_int, c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
        POINTER(c_int), POINTER(c_int),
        POINTER(c_float), POINTER(c_int)
    ]

    n_len, m_len = X.shape

    h_arr = np.arange(1, k * m_len + 1, dtype=np.float32).reshape((k, m_len))
    w_arr = np.arange(1, n_len * k + 1, dtype=np.float32).reshape((n_len, k))

    fn(X, c_int(n_len), c_int(m_len), c_int(k),
       w_arr, h_arr,
       pointer(c_int(max_iter)), pointer(c_int(tp)),
       pointer(c_float(thresh)), pointer(c_int(norm)))

    return h_arr, w_arr
