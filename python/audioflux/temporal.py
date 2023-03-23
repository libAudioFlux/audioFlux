import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel

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
    >>> energy_arr, rms_arr, zero_cross_arr, m_arr = temp_obj.get_data(audio_arr)
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
        self._is_created = True

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

    def get_data(self, data_arr):
        """
           Get energy/rms/zeroCrossRate feature

           Parameters
           ----------
            data_arr: np.ndarray [shape=(..., n)]
                Input audio data

           Returns
           -------
           energy_arr: np.ndarray [shape=(..., time)]
               energy feature

           rms_arr: np.ndarray [shape=(..., time)]
               rms feature

           zcr_arr: np.ndarray [shape=(..., time)]
               zero Cross Rate feature

           m_arr: np.ndarray [shape=(..., time, frame)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        data_length = data_arr.shape[-1]
        if data_arr.ndim == 1:
            self._temporal(data_arr)
            energy_arr, rms_arr, zcr_arr, m_arr = self._get_data(data_length)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            energy_arr = []
            rms_arr = []
            zcr_arr = []
            m_arr = []

            for i in range(channel_num):
                self._temporal(data_arr[i])
                _energy_arr, _rms_arr, _zcr_arr, _m_arr = self._get_data(data_length)
                energy_arr.append(_energy_arr)
                rms_arr.append(_rms_arr)
                zcr_arr.append(_zcr_arr)
                m_arr.append(_m_arr)

            energy_arr = np.stack(energy_arr, axis=0)
            rms_arr = np.stack(rms_arr, axis=0)
            zcr_arr = np.stack(zcr_arr, axis=0)
            m_arr = np.stack(m_arr, axis=0)

            energy_arr = revoke_channel(energy_arr, o_channel_shape, 1)
            rms_arr = revoke_channel(rms_arr, o_channel_shape, 1)
            zcr_arr = revoke_channel(zcr_arr, o_channel_shape, 1)
            m_arr = revoke_channel(m_arr, o_channel_shape, 2)
        return energy_arr, rms_arr, zcr_arr, m_arr

    def _temporal(self, data_arr):
        fn = self._lib['temporalObj_temporal']
        fn.argtypes = [POINTER(OpaqueTemporal),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       ]
        data_length = data_arr.size
        fn(self._obj, data_arr, c_int(data_length))

    def _get_data(self, data_length):
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
                         dtype=np.float32).reshape(time_length, self.frame_length)

        return energy_arr, rms_arr, zero_cross_arr, m_arr

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['temporal_free']
            free_fn.argtypes = [POINTER(OpaqueTemporal)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
