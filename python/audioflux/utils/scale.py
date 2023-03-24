from ctypes import c_int
import numpy as np

from audioflux.fftlib import get_fft_lib

__all__ = [
    'min_max_scale',
    'stand_scale',
    'max_abs_scale',
    'robust_scale',
    'center_scale',
    'mean_scale',
    'arctan_scale',
]


def min_max_scale(X):
    """
    min max scale
    
    Parameters
    ----------
    X: np.ndarray [shape=(n_samples, n_features)]
        Input array.

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.
        
    """""
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_minMaxScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def stand_scale(X, tp: int = 1):
    """
    stand scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array.

    tp: int
        - 0 sample variance
        - 1 population variance

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_standScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), c_int(tp), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def max_abs_scale(X):
    """
    max abs scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array.

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_maxAbsScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def robust_scale(X):
    """
    robust scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_robustScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def center_scale(X):
    """
    center scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.

    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_centerScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def mean_scale(X):
    """
    mean scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.

    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_meanScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)


def arctan_scale(X):
    """
    arctan scale

    Parameters
    ----------
    X: np.ndarray (n_samples, n_features)
        Input array

    Returns
    -------
    out: np.ndarray [shape=(n_samples, n_features)]
        Transformed array.

    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_arctanScale']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    sample_len, feature_len = X.shape
    ret_arr = np.zeros_like(X, dtype=X.dtype)
    for i in range(feature_len):
        feature_arr = X[:, i]
        feature_arr = np.ascontiguousarray(feature_arr)
        feature_scale_arr = np.zeros_like(feature_arr, dtype=feature_arr.dtype)
        fn(feature_arr, c_int(sample_len), feature_scale_arr)
        ret_arr[:, i] = feature_scale_arr
    return np.ascontiguousarray(ret_arr)
