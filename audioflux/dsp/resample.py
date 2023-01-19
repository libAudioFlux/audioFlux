import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.base import Base
from audioflux.type import ResampleQualityType, WindowType
from audioflux.utils import check_audio

__all__ = ["Resample", "WindowResample"]


class OpaqueResample(Structure):
    _fields_ = []


class ResampleBase(Base):
    def __init__(self):
        super(ResampleBase, self).__init__(pointer(OpaqueResample()))

    def set_samplate(self, source_rate, target_rate):
        """
        Set samplate

        Parameters
        ----------
        source_rate: int
        target_rate: int
        """
        fn = self._lib['resampleObj_setSamplate']
        fn.argtypes = [POINTER(OpaqueResample), c_int, c_int]
        fn(self._obj, c_int(source_rate), c_int(target_rate))

    def enable_continue(self, flag):
        """
        Enable continue

        Parameters
        ----------
        flag: int
        """
        fn = self._lib['resampleObj_enableContinue']
        fn.argtypes = [POINTER(OpaqueResample), c_int]
        fn(self._obj, c_int(flag))

    def cal_data_length(self, data_length):
        """
        Compute the data length

        Parameters
        ----------
        data_length: int

        Returns
        -------
        out: int
        """
        fn = self._lib['resampleObj_calDataLength']
        fn.argtypes = [POINTER(OpaqueResample), c_int]
        fn.restype = c_int
        fn(self._obj, c_int(data_length))

    def debug(self):
        fn = self._lib['resampleObj_debug']
        fn.argtypes = [POINTER(OpaqueResample)]
        fn(self._obj)

    def resample(self, data_arr):
        """
        Resample for audio data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]
            Input audio array

        Returns
        -------
        out: np.ndarray [shape=(n,)]
            Audio data after resampling
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)

        fn = self._lib['resampleObj_resample']
        fn.argtypes = [POINTER(OpaqueResample),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')]
        data_len = data_arr.shape[0]

        ret_arr = np.zeros(data_len * 5, dtype=np.float32)

        new_arr_len = fn(self._obj, data_arr, data_len, ret_arr)
        return ret_arr[:new_arr_len]

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['resampleObj_free']
            free_fn.argtypes = [POINTER(OpaqueResample)]
            free_fn.restype = c_void_p
            free_fn(self._obj)


class Resample(ResampleBase):
    """
    Calculate resampling

    Parameters
    ----------
    qual_type: ResampleQualityType
    is_scale: bool
        Whether to use scale

    Examples
    --------
    >>> import audioflux as af
    >>> from audioflux.type import ResampleQualityType
    >>> audio_data, sr = af.read(af.utils.sample_path('220'))
    >>>
    >>> resample_obj = af.Resample(qual_type=ResampleQualityType.BEST, is_scale=False)
    >>> new_audio_data = resample_obj.resample(audio_data)
    >>> new_audio_data
    array([ 7.7589704e-07, -1.4843295e-06,  2.4068654e-06, ...,
            2.6401449e-03,  3.1324192e-03,  3.3868181e-03], dtype=float32)
    """

    def __init__(self, qual_type=ResampleQualityType.BEST, is_scale=False):
        super(Resample, self).__init__()

        self.qual_type = qual_type
        self.is_scale = is_scale
        self.is_continue = False

        fn = self._lib['resampleObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueResample)),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           None if self.qual_type is None else pointer(c_int(self.qual_type.value)),
           pointer(c_int(int(self.is_scale))),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class WindowResample(ResampleBase):
    """
    Calculate window resample

    Parameters
    ----------
    zero_num: int
    nbit: int
    win_type: WindowType
    value: float
    roll_off: float
    is_scale: bool
        Whether to use scale

    Examples
    --------
    >>> import audioflux as af
    >>> from audioflux.type import ResampleQualityType
    >>> audio_data, sr = af.read(af.utils.sample_path('220'))
    >>>
    >>> resample_obj = af.WindowResample()
    >>> new_audio_data = resample_obj.resample(audio_data)
    >>> new_audio_data
    array([ 1.0972841e-06, -2.0991588e-06,  3.4038217e-06, ...,
            3.7337288e-03,  4.4299099e-03,  4.7896840e-03], dtype=float32)
    """

    def __init__(self, zero_num=64, nbit=9, win_type=WindowType.HANN,
                 value=None, roll_off=0.945, is_scale=False):
        super(WindowResample, self).__init__()

        self.zero_num = zero_num
        self.nbit = nbit
        self.win_type = win_type
        self.value = value
        self.roll_off = roll_off
        self.is_scale = is_scale
        self.is_continue = False

        fn = self._lib['resampleObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueResample)),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_float),
                       POINTER(c_float),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           pointer(c_int(self.zero_num)),
           pointer(c_int(self.nbit)),
           pointer(c_int(self.win_type.value)),
           None if self.value is None else pointer(c_float(self.value)),
           pointer(c_float(self.roll_off)),
           pointer(c_int(int(self.is_scale))),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True
