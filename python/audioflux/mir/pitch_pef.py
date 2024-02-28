import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel

__all__ = ["PitchPEF"]


class OpaquePitchPEF(Structure):
    _fields_ = []


class PitchPEF(Base):
    """
    Pitch Estimation Filter(PEF). A pitch estimation filter is designed, and pitch is estimated by performing
    cross-correlation operations in the frequency domain. [#]_

    .. [#] Gonzalez, Sira, and Mike Brookes. "A Pitch Estimation Filter robust to high levels of
           noise (PEFAC)." 19th European Signal Processing Conference. Barcelona, 2011, pp. 451–455.

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. Default is `32.0`.

    high_fre: float
        Highest frequency. Default is `2000.0`.

    cut_fre: float
        Cut frequency. Default is `4000.0`, and must be greater than `high_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    slide_length: int
        Window sliding length.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    alpha: float, > 0
        alpha. Default if `10.0`..

    beta: float, 0~1
        beta. Default if `0.5`..

    gamma: float, > 1
        gamma. Default if `1.8`.

    Examples
    --------

    Read voice audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('voice')
    >>> audio_arr, sr = af.read(audio_path)

    Extract pitch

    >>> pitch_obj = af.PitchPEF(samplate=sr)
    >>> fre_arr = pitch_obj.pitch(audio_arr)

    Show pitch plot

    >>> import matplotlib.pyplot as plt
    >>> times = np.arange(fre_arr.shape[-1]) * (pitch_obj.slide_length / sr)
    >>> fig, ax = plt.subplots(nrows=2, sharex=True)
    >>> ax[0].set_title('PitchPEF')
    >>> af.display.fill_wave(audio_arr, samplate=sr, axes=ax[0])
    >>> ax[1].scatter(times, fre_arr, s=2)
    """

    def __init__(self, samplate=32000, low_fre=32.0, high_fre=2000.0, cut_fre=4000.0,
                 radix2_exp=12, slide_length=1024, window_type=WindowType.HAMM,
                 alpha=10.0, beta=0.5, gamma=1.8):
        super(PitchPEF, self).__init__(pointer(OpaquePitchPEF()))

        if low_fre >= high_fre:
            raise ValueError(f'`low_fre` must be smaller than `high_fre`')
        if high_fre >= cut_fre:
            raise ValueError(f'`high_fre` must be smaller than `cut_fre`')
        if alpha <= 0:
            raise ValueError(f'`alpha` must be greater than 0.')
        if beta < 0 or beta > 1:
            raise ValueError(f'`beta` must be between 0 and 1.')
        if gamma <= 1:
            raise ValueError(f'`gamma` must be greater than 1.')

        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.cut_fre = cut_fre
        self.radix2_exp = radix2_exp
        self.slide_length = slide_length
        self.window_type = window_type
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.is_continue = False

        fn = self._lib['pitchPEFObj_new']
        fn.argtypes = [POINTER(POINTER(OpaquePitchPEF)), POINTER(c_int),
                       POINTER(c_float), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int),
                       POINTER(c_float), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int)]

        fn(self._obj,
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_float(self.cut_fre)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.window_type.value)),
           pointer(c_float(self.alpha)),
           pointer(c_float(self.beta)),
           pointer(c_float(self.gamma)),
           pointer(c_int(int(self.is_continue))))
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
        fn = self._lib['pitchPEFObj_calTimeLength']
        fn.argtypes = [POINTER(OpaquePitchPEF), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def set_filter_params(self, alpha, beta, gamma):
        """
        Set filter params

        Parameters
        ----------
        alpha: float
            alpha

        beta: float
            beta

        gamma: float
            gamma

        """
        if alpha <= 0:
            raise ValueError(f'`alpha` must be greater than 0.')
        if beta < 0 or beta > 1:
            raise ValueError(f'`beta` must be between 0 and 1.')
        if gamma <= 1:
            raise ValueError(f'`gamma` must be greater than 1.')

        fn = self._lib['pitchPEFObj_setFilterParams']
        fn.argtypes = [POINTER(OpaquePitchPEF), c_float, c_float, c_float]
        fn(self._obj, c_float(alpha), c_float(beta), c_float(gamma))

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

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
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['pitchPEFObj_pitch']
        fn.argtypes = [POINTER(OpaquePitchPEF),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        data_len = data_arr.shape[-1]
        time_length = self.cal_time_length(data_len)

        if data_arr.ndim == 1:
            fre_arr = np.zeros(time_length, dtype=np.float32)
            fn(self._obj, data_arr, c_int(data_len), fre_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            size = (channel_num, time_length)
            fre_arr = np.zeros(size, dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], c_int(data_len), fre_arr[i])

            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)
        return fre_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['pitchPEFObj_free']
            fn.argtypes = [POINTER(OpaquePitchPEF)]
            fn.restype = c_void_p
            fn(self._obj)
