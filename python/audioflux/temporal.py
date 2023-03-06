import ctypes

import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType

__all__ = ['Temporal']


class OpaqueTemporal(Structure):
    _fields_ = []


class Temporal(Base):
    """
    Temporal feature

    Parameters
    ----------
    frame_length: int
        frame length

    slide_length: int
        sliding length

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Create Temporal and extract feature

    >>> temp_obj = af.Temporal(frame_length=2048, slide_length=512)
    >>> temp_obj.temporal(audio_arr)
    >>> energy_arr, rms_arr, zero_cross_arr, m_arr = temp_obj.get_data(len(audio_arr))
    """

    def __init__(self, frame_length=2048, slide_length=512, window_type=WindowType.HANN):
        super(Temporal, self).__init__(pointer(OpaqueTemporal()))

        self.frame_length = frame_length
        self.slide_length = slide_length
        self.window_type = window_type

        fn = self._lib['temporalObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueTemporal)),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int)]
        fn(self._obj,
           pointer(c_int(self.frame_length)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.window_type.value)))
        self._is_created = False

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['temporalObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueTemporal), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def temporal(self, data_arr):
        """
        set audio data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]
            Input audio data
        """
        data_arr = data_arr.astype(dtype=np.float32)

        fn = self._lib['temporalObj_temporal']
        fn.argtypes = [POINTER(OpaqueTemporal),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       ]
        data_length = data_arr.size
        fn(self._obj, data_arr, c_int(data_length))

    def get_data(self, data_length) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Get energy/rms/zeroCrossRate feature

        Parameters
        ----------
        data_length: int
            The length of the data passed in by the temporal method

        Returns
        -------
        energy_arr: np.ndarray [shape=(time,)]
            energy feature

        rms_arr: np.ndarray [shape=(time,)]
            rms feature

        zero_cross_arr: np.ndarray [shape=(time,)]
            zero Cross Rate feature

        m_arr: np.ndarray [shape=(...,time)]
        """
        fn = self._lib['temporalObj_getData']
        fn.argtypes = [POINTER(OpaqueTemporal),
                       POINTER(POINTER(c_float)),
                       POINTER(POINTER(c_float)),
                       POINTER(POINTER(c_float)),
                       POINTER(POINTER(c_float)),
                       ]
        time_length = self.cal_time_length(data_length)

        pp_energy_arr = pointer(pointer(c_float()))
        pp_rms_arr = pointer(pointer(c_float()))
        pp_zero_cross_arr = pointer(pointer(c_float()))
        pp_m_arr = pointer(pointer(c_float()))

        fn(self._obj,
           pp_energy_arr,
           pp_rms_arr,
           pp_zero_cross_arr,
           pp_m_arr,
           )

        energy_arr = np.array([pp_energy_arr.contents[x] for x in range(time_length)], dtype=np.float32)
        rms_arr = np.array([pp_rms_arr.contents[x] for x in range(time_length)], dtype=np.float32)
        zero_cross_arr = np.array([pp_zero_cross_arr.contents[x] for x in range(time_length)], dtype=np.float32)
        m_arr = np.array([pp_m_arr.contents[x] for x in range(time_length * self.frame_length)],
                         dtype=np.float32).reshape(
            time_length, self.frame_length)

        return energy_arr, rms_arr, zero_cross_arr, m_arr

    def ezr(self, data_length, gamma):
        """
        Get ezr feature

        Parameters
        ----------
        data_length: int
            The length of the data passed in by the temporal method

        gamma: float
            gamma value

        Returns
        -------
        out: np.ndarray [shape=(time,)]
            ezr feature
        """

        fn = self._lib['temporalObj_ezr']
        fn.argtypes = [POINTER(OpaqueTemporal),
                       c_float,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]
        time_len = self.cal_time_length(data_length)
        ret_arr = np.zeros(time_len, dtype=np.float32)
        fn(self._obj, c_float(gamma), ret_arr)
        return ret_arr

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['temporal_free']
            free_fn.argtypes = [POINTER(OpaqueTemporal)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
