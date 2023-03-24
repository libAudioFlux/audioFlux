import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import PitchType
from audioflux.utils import check_audio, format_channel, revoke_channel, note_to_hz

__all__ = ["Pitch"]


class OpaquePitch(Structure):
    _fields_ = []


class Pitch(Base):
    """
    Pitch - YIN, STFT, etc algorithm

    Parameters
    ----------
    pitch_type: PitchType
        Pitch type

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    high_fre: float
        Highest frequency.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    slide_length: int
        Window sliding length.

    auto_length: int
        Auto length

    Examples
    --------

    Get a 220Hz's audio file

    >>> import audioflux as af
    >>> audio_arr, sr = af.read(af.utils.sample_path('220'))
    # >>> audio_arr = audio_arr[:8192]

    Create Pitch object and get frequency

    >>> from audioflux.type import PitchType
    >>> obj = af.Pitch(pitch_type=PitchType.YIN)
    >>> fre_arr, value_arr1, value_arr2 = obj.pitch(audio_arr)

    Display plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_wave, fill_plot
    >>> import numpy as np
    >>> audio_len = audio_arr.shape[-1]
    >>> fig, axes = plt.subplots(nrows=2)
    >>> fill_wave(audio_arr, samplate=sr, axes=axes[0])
    >>>
    >>> ax = fill_plot(np.arange(len(fre_arr)), fre_arr, label='fre', axes=axes[1])
    >>> ax.set_ylabel('frequency(Hz)')
    """

    def __init__(self, pitch_type=None, samplate=32000,
                 low_fre=note_to_hz('A0'), high_fre=note_to_hz('C7'),
                 radix2_exp=12, slide_length=1024, auto_length=2048):
        super(Pitch, self).__init__(pointer(OpaquePitch()))

        self.pitch_type = pitch_type
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.radix2_exp = radix2_exp
        self.slide_length = slide_length
        self.auto_length = auto_length
        self.is_continue = False

        fn = self._lib['pitchObj_new']
        fn.argtypes = [POINTER(POINTER(OpaquePitch)),
                       POINTER(c_int), POINTER(c_int),
                       POINTER(c_float), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int),
                       POINTER(c_int), POINTER(c_int)]

        fn(self._obj,
           None if self.pitch_type is None else pointer(c_int(self.pitch_type.value)),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.auto_length)),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True

    def set_thresh(self, thresh):
        """
        Set thresh

        Parameters
        ----------
        thresh: float
        """
        fn = self._lib['pitchObj_setThresh']
        fn.argtypes = [POINTER(OpaquePitch), c_float]
        fn(self._obj, c_float(thresh))

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        - ``fft_length = 2 ** radix2_exp``
        - ``(data_length - fft_length) / slide_length + 1``

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['pitchObj_calTimeLength']
        fn.argtypes = [POINTER(OpaquePitch), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def pitch(self, data_arr):
        """
        Compute pitch

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio array

        Returns
        -------
        fre_arr: np.ndarray [shape=(..., time)]
        value1_arr: np.ndarray [shape=(..., time)]
        value2_arr: np.ndarray [shape=(..., time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['pitchObj_pitch']
        fn.argtypes = [POINTER(OpaquePitch),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        data_len = data_arr.shape[-1]
        time_length = self.cal_time_length(data_len)

        if data_arr.ndim == 1:
            fre_arr = np.zeros(time_length, dtype=np.float32)
            value1_arr = np.zeros(time_length, dtype=np.float32)
            value2_arr = np.zeros(time_length, dtype=np.float32)
            fn(self._obj, data_arr, c_int(data_len), fre_arr, value1_arr, value2_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            size = (channel_num, time_length)
            fre_arr = np.zeros(size, dtype=np.float32)
            value1_arr = np.zeros(size, dtype=np.float32)
            value2_arr = np.zeros(size, dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], c_int(data_len), fre_arr[i], value1_arr[i], value2_arr[i])

            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)
            value1_arr = revoke_channel(value1_arr, o_channel_shape, 1)
            value2_arr = revoke_channel(value2_arr, o_channel_shape, 1)
        return fre_arr, value1_arr, value2_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['pitchObj_free']
            fn.argtypes = [POINTER(OpaquePitch)]
            fn.restype = c_void_p
            fn(self._obj)
