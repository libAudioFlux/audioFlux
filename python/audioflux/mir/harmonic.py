import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.type import WindowType
from audioflux.base import Base
from audioflux.utils import check_audio, format_channel, revoke_channel

__all__ = ["Harmonic"]


class OpaqueHarmonic(Structure):
    _fields_ = []


class Harmonic(Base):
    """
    Harmonics

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    slide_length: int
        Window sliding length.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    low_fre: float
        Lowest frequency. Default is `27.0`.

    high_fre: float
        Highest frequency. Default is `4000.0`.

    Examples
    --------

    Read guitar_chord2 audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('guitar_chord2')
    >>> audio_arr, sr = af.read(audio_path)

    Compute harmonic_count

    >>> hr_obj = af.Harmonic(radix2_exp=12, samplate=sr, slide_length=1024, window_type=af.type.WindowType.HAMM)
    >>> count_arr = hr_obj.harmonic_count(audio_arr, 82, 900)

    Show harmonic_count plot

    >>> import matplotlib.pyplot as plt
    >>> fig, axes = plt.subplots(2, figsize=(16, 9), sharex=True)
    >>> af.display.fill_wave(audio_arr, axes=axes[0])
    >>> times = np.arange(count_arr.shape[-1]) * (hr_obj.slide_length / sr)
    >>> af.display.fill_plot(times, count_arr, axes=axes[1], label='Harmonic Count')
    """

    def __init__(self, radix2_exp=12, samplate=32000, slide_length=1024,
                 window_type=WindowType.HAMM, low_fre=27., high_fre=4000.):
        super(Harmonic, self).__init__(pointer(OpaqueHarmonic()))

        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.slide_length = slide_length
        self.window_type = window_type
        self.low_fre = low_fre
        self.high_fre = high_fre

        self.fft_length = 1 << self.radix2_exp

        fn = self._lib['harmonicObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueHarmonic)),
                       POINTER(c_int), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int)]
        fn(self._obj,
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)))
        self._is_created = True

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        - ``fft_length = 2 ** radix2_exp``
        - ``(data_length - fft_length) // slide_length + 1``

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['harmonicObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueHarmonic), c_int]
        return fn(self._obj, c_int(data_length))

    def _exec(self, data_arr):
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        exec_fn = self._lib['harmonicObj_exec']
        exec_fn.argtypes = [
            POINTER(OpaqueHarmonic),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int
        ]
        exec_fn(self._obj, data_arr, c_int(data_arr.shape[-1]))

    def _harmonic_count(self, data_length, low_fre, high_fre):
        hc_fn = self._lib['harmonicObj_harmonicCount']
        hc_fn.argtypes = [
            POINTER(OpaqueHarmonic),
            c_float,
            c_float,
            np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_time = self.cal_time_length(data_length)
        count_arr = np.zeros(n_time, dtype=np.int32)
        hc_fn(self._obj, c_float(low_fre), c_float(high_fre), count_arr)
        return count_arr

    def harmonic_count(self, data_arr, low_fre, high_fre):
        """

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input data array.

        low_fre: float
            Lowest frequency.
        high_fre: float
            Highest frequency.

        Returns
        -------
        out: np.ndarray [shape=(..., time)]
            Harmonic count array.
        """

        if self.low_fre > low_fre:
            raise ValueError(f'low_fre must be greater than or equal to {self.low_fre}')
        if self.high_fre < high_fre:
            raise ValueError(f'high_fre must be less than or equal to {self.high_fre}')
        if low_fre > high_fre:
            raise ValueError(f'low_fre must be less than or equal to high_fre')

        data_len = data_arr.shape[-1]
        if data_arr.ndim == 1:
            self._exec(data_arr)
            count_arr = self._harmonic_count(data_len, low_fre, high_fre)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]
            count_arr = []
            for i in range(channel_num):
                self._exec(data_arr[i])
                _count_arr = self._harmonic_count(data_len, low_fre, high_fre)
                count_arr.append(_count_arr)
            count_arr = np.stack(count_arr, axis=0)
            count_arr = revoke_channel(count_arr, o_channel_shape, 1)
        return count_arr

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['harmonicObj_free']
            free_fn.argtypes = [POINTER(OpaqueHarmonic)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
