import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel, note_to_hz

__all__ = ["HarmonicRatio"]


class OpaqueHarmonicRatio(Structure):
    _fields_ = []


class HarmonicRatio(Base):
    """
    Harmonic Ratio

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    """

    def __init__(self, samplate=32000, low_fre=note_to_hz('C1'), radix2_exp=12,
                 window_type=WindowType.HAMM, slide_length=1024):
        super(HarmonicRatio, self).__init__(pointer(OpaqueHarmonicRatio()))

        self.samplate = samplate
        self.low_fre = low_fre
        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length

        fn = self._lib['harmonicRatioObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueHarmonicRatio)),
                       POINTER(c_int), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int),
                       POINTER(c_int)]

        fn(self._obj,
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)))
        self._is_created = True

    def cal_time_length(self, data_length):
        """
        Compute the time length

        Parameters
        ----------
        data_length: int
            Input array length

        Returns
        -------
        out: int
        """
        fn = self._lib['harmonicRatioObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueHarmonicRatio), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def harmonic_ratio(self, data_arr):
        """
        Compute harmonic ratio

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data array.

        Returns
        -------
        out: np.ndarray [shape=(..., time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['harmonicRatioObj_harmonicRatio']
        fn.argtypes = [POINTER(OpaqueHarmonicRatio),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        data_len = data_arr.shape[-1]
        time_length = self.cal_time_length(data_len)

        if data_arr.ndim == 1:
            ret_arr = np.zeros(time_length, dtype=np.float32)
            fn(self._obj, data_arr, c_int(data_len), ret_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            ret_arr = np.zeros((channel_num, time_length), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], c_int(data_len), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['harmonicRatioObj_free']
            fn.argtypes = [POINTER(OpaqueHarmonicRatio)]
            fn.restype = c_void_p
            fn(self._obj)
