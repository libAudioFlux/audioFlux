import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import XcorrNormalType

__all__ = ['Xcorr']


class OpaqueXcorr(Structure):
    _fields_ = []


class Xcorr(Base):
    """
    Cross-correlation of two 1-dimensional sequences.

    Examples
    --------
    >>> import audioflux as af
    >>> import numpy as np
    >>>
    >>> arr1 = np.array([1, 2, 3], dtype=np.float32)
    >>> arr2 = np.array([0, 1, 0.5], dtype=np.float32)
    >>> xcorr = af.Xcorr()
    >>> xcorr_arr, max_val = xcorr.xcorr(arr1, arr2)
    >>> print(xcorr_arr)
    >>> print(max_val)
    """

    def __init__(self):
        super(Xcorr, self).__init__(pointer(OpaqueXcorr()))

        fn = self._lib['xcorrObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueXcorr))]
        fn(self._obj)
        self._is_created = True

    def xcorr(self, data_arr1, data_arr2=None, xcorr_normal_type=XcorrNormalType.NONE):
        """
        Calculate cross-correlation (xcorr).

        Parameters
        ----------
        data_arr1: np.ndarray [shape=(n)]
            data array 1

        data_arr2: np.ndarray [shape=(n)] or None
            data array 2

            If it is None, calculate auto-correlation.

        xcorr_normal_type: XcorrNormalType
            Normalization type

            See: `type.XcorrNormalType`

        Returns
        -------
        arr: np.ndarray [shape=(2 * n - 1)]
            X-correlation array
        max_val: float
            Maximum value
        """
        data_arr1 = np.asarray(data_arr1, dtype=np.float32, order='C')
        data_arr2 = None if data_arr2 is None else np.asarray(data_arr2, dtype=np.float32, order='C')

        if data_arr1.ndim != 1:
            raise ValueError(f"data_arr1[ndim={data_arr1.ndim}] must be a 1D array")
        if data_arr2 is not None and data_arr2.ndim != 1:
            raise ValueError(f"data_arr2[ndim={data_arr2.ndim}] must be a 1D array")
        if data_arr2 is not None and data_arr1.shape != data_arr2.shape:
            raise ValueError(f"data_arr1.shape={data_arr1.shape} must be equal to data_arr2.shape={data_arr2.shape}")

        data_length = data_arr1.shape[0]

        fn = self._lib['xcorrObj_xcorr']
        fn.argtypes = [
            POINTER(OpaqueXcorr),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            POINTER(c_int) if data_arr2 is None else np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,
                                                                            flags='C_CONTIGUOUS'),
            c_int, POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            POINTER(c_float)
        ]

        arr = np.zeros(2 * data_length - 1, dtype=np.float32)
        max_val = pointer(c_float())
        fn(self._obj, data_arr1, data_arr2, c_int(data_length),
           pointer(c_int(xcorr_normal_type.value)), arr, max_val)

        return arr, max_val.contents.value

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['xcorrObj_free']
            free_fn.argtypes = [POINTER(OpaqueXcorr)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
