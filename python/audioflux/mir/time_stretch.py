import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel

__all__ = ["TimeStretch"]


class OpaqueTimeStretch(Structure):
    _fields_ = []


class TimeStretch(Base):
    """
    Time stretch algorithm

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    slide_length: int
        Window sliding length.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Examples
    --------

    Read voice audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('voice')
    >>> audio_arr, sr = af.read(audio_path)

    Compute the time stretch

    >>> time_stretch_obj = af.TimeStretch(radix2_exp=12, window_type=af.type.WindowType.HANN, slide_length=1024)
    >>> new_audio_arr = time_stretch_obj.time_stretch(audio_arr, 0.5)
    >>> # af.write('./audio_arr_0_5.wav', new_audio_arr, sr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(nrows=2, sharex=True)
    >>> ax = af.display.fill_wave(audio_arr, sr, axes=axes[0])
    >>> ax.set_title('Original')
    >>> ax = af.display.fill_wave(new_audio_arr, sr, axes=axes[1])
    >>> ax.set_title('0.5x TimeStretch')

    """

    def __init__(self, radix2_exp=12, slide_length=1024, window_type=WindowType.HANN):
        super(TimeStretch, self).__init__(pointer(OpaqueTimeStretch()))

        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length

        fn = self._lib['timeStretchObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueTimeStretch)),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]

        fn(self._obj,
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.window_type.value)))
        self._is_created = True

    def cal_data_capacity(self, rate, data_length):
        """
        Calculate the data capacity.

        Parameters
        ----------
        rate: float
            Time stretch rate

        data_length: int
            Input array length

        Returns
        -------
        out: int
        """
        fn = self._lib['timeStretchObj_calDataCapacity']
        fn.argtypes = [POINTER(OpaqueTimeStretch), c_float, c_int]
        fn.restype = c_int
        return fn(self._obj, c_float(rate), c_int(data_length))

    def time_stretch(self, data_arr, rate):
        """
        Compute the time stretch

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Audio data array

        rate: float
            Time stretch rate

        Returns
        -------
        arr: np.ndarray [shape=(..., n1)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['timeStretchObj_timeStretch']
        fn.argtypes = [
            POINTER(OpaqueTimeStretch),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]

        data_len = data_arr.shape[-1]
        new_data_len = self.cal_data_capacity(rate, data_len)
        rate_c = c_float(rate)
        data_len_c = c_int(data_len)

        if data_arr.ndim == 1:
            new_arr = np.zeros(new_data_len, dtype=np.float32)
            fn(self._obj, rate_c, data_arr, data_len_c, new_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            new_arr = np.zeros((channel_num, new_data_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, rate_c, data_arr[i], data_len_c, new_arr[i])

            new_arr = revoke_channel(new_arr, o_channel_shape, 1)

        return new_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['timeStretchObj_free']
            fn.argtypes = [POINTER(OpaqueTimeStretch)]
            fn.restype = c_void_p
            fn(self._obj)
