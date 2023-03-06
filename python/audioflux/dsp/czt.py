import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.base import Base
from audioflux.utils import check_audio

__all__ = ["CZT"]


class OpaqueCZT(Structure):
    _fields_ = []


class CZT(Base):
    """
    Chirp Z-transform (CZT)

    Parameters
    ----------
    radix2_exp: int
        ``fft_length = 2 ** radix2_exp``

    Examples
    --------
    >>> import numpy as np
    >>> from audioflux import CZT
    >>> n = np.arange(4096)
    >>> c_obj = CZT(3)
    >>> czt_arr = c_obj.czt(n, 0.15, 0.25)
    >>> czt_arr = np.abs(czt_arr)
    >>> czt_arr
    array([7.4057508, 6.685393 , 6.4312935, ..., 0.       , 0.       ,
           0.       ], dtype=float32)
    """

    def __init__(self, radix2_exp):
        super(CZT, self).__init__(pointer(OpaqueCZT()))

        self.radix2_exp = radix2_exp

        fn = self._lib['cztObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueCZT)), c_int]
        fn(self._obj, c_int(self.radix2_exp))
        self._is_created = True

    def czt(self, data_arr, low_w, high_w):
        """
        Compute the czt

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]
            Input audio 1D array
        low_w: float
        high_w: float

        Returns
        -------
        out: np.ndarray [shape=(n_time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)

        czt_fn = self._lib['cztObj_czt']
        czt_fn.argtypes = [POINTER(OpaqueCZT),
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                           c_float, c_float,
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')]

        real_arr = data_arr.real
        imag_arr = data_arr.imag

        data_length = data_arr.size
        size = data_length * 2

        ret_real_arr = np.zeros(size, dtype=np.float32)
        ret_imag_arr = np.zeros(size, dtype=np.float32)

        czt_fn(self._obj, real_arr, imag_arr,
               c_float(low_w), c_float(high_w),
               ret_real_arr, ret_imag_arr)
        return ret_real_arr + (ret_imag_arr * 1j)

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['cztObj_free']
            free_fn.argtypes = [POINTER(OpaqueCZT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
