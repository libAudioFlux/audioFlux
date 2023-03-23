import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.type import WaveletContinueType, SpectralFilterBankScaleType, get_wavelet_default_gamma_beta
from audioflux.base import Base
from audioflux.utils import check_audio, check_audio_length, format_channel, revoke_channel, note_to_hz

__all__ = ['WSST']


class OpaqueWSST(Structure):
    _fields_ = []


class WSST(Base):
    """
    Wavelet synchrosqueezed transform (WSST)

    Parameters
    ----------
    num: int
        Number of frequency bins to generate

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio

    low_fre: float
        Lowest frequency

        - Linear/Linsapce/Mel/Bark/Erb, low_fre>=0. `default: 0.0`
        - Octave/Log, low_fre>=32.703. `default: 32.703(C1)`

    high_fre: float
        Highest frequency.

        Default is `16000(samplate/2)`
        Octave is not provided, it is based on musical pitch

    bin_per_octave: int
        Number of bins per octave.

    scale_type: SpectralFilterBankScaleType
        Spectral filter bank type. It determines the type of spectrogram.

        See: `type.SpectralFilterBankScaleType`

    wavelet_type: WaveletContinueType
        Wavelet continue type

        See: `type.WaveletContinueType`

        .. note::
            Default gamma/beta values for different wavelet_types:

            - morse: gamma=3 beta=20
            - morlet: gamma=6 beta=2
            - bump: gamma=5 beta=0.6
            - paul: gamma 4
            - dog: gamma 2 beta 2; must even
            - mexican: beta 2
            - hermit: gamma 5 beta 2
            - ricker: gamma 4

    gamma: float or None
        gamma value

    beta: float or None
        beta value

    thresh: float
        thresh

    is_padding: bool
        Whether to use padding.

    See Also
    --------
    Reassign
    Synsq

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)
    >>> # WSST can only input fft_length data
    >>> # For radix2_exp=12, then fft_length=4096
    >>> audio_arr = audio_arr[..., :4096]

    Create WSST object of mel

    >>> from audioflux.type import SpectralFilterBankScaleType, WaveletContinueType
    >>> from audioflux.utils import note_to_hz
    >>> obj = af.WSST(num=128, radix2_exp=12, samplate=sr,
    >>>               bin_per_octave=12, wavelet_type=WaveletContinueType.MORSE,
    >>>               scale_type=SpectralFilterBankScaleType.MEL, is_padding=False)

    Extract spectrogram

    >>> import numpy as np
    >>> wsst_spec_arr, cwt_spec_arr = obj.wsst(audio_arr)
    >>> wsst_spec_arr = np.abs(wsst_spec_arr)
    >>> cwt_spec_arr = np.abs(cwt_spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> # Show CWT
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cwt_spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='CWT-Mel')
    >>> fig.colorbar(img, ax=ax)
    >>>
    >>> # Show WSST
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(wsst_spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='WSST-Mel')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, num=84, radix2_exp=12, samplate=32000, low_fre=None, high_fre=None,
                 bin_per_octave=12, wavelet_type=WaveletContinueType.MORLET,
                 scale_type=SpectralFilterBankScaleType.OCTAVE,
                 gamma=None, beta=None, thresh=0.001, is_padding=True):

        super(WSST, self).__init__(pointer(OpaqueWSST()))

        # check BPO
        if scale_type == SpectralFilterBankScaleType.OCTAVE and bin_per_octave < 1:
            raise ValueError(f'bin_per_octave={bin_per_octave} must be a positive integer')

        if low_fre is None:
            if scale_type in (SpectralFilterBankScaleType.OCTAVE,
                              SpectralFilterBankScaleType.LOG):
                low_fre = note_to_hz('C1')  # 32.703
            else:
                low_fre = 0.0

        if high_fre is None:
            high_fre = samplate / 2

        # check low_fre
        if scale_type in (SpectralFilterBankScaleType.OCTAVE,
                          SpectralFilterBankScaleType.LOG) \
                and low_fre < round(note_to_hz('C1'), 3):
            # Octave/Log >= 32.703
            raise ValueError(f'{scale_type.name} low_fre={low_fre} must be greater than or equal to 32.703')
        if low_fre < 0:
            # linear/linspace/mel/bark/erb low_fre>=0
            raise ValueError(f'{scale_type.name} low_fre={low_fre} must be a non-negative number')

        default_gamma, default_beta = get_wavelet_default_gamma_beta(wavelet_type)
        if gamma is None:
            gamma = default_gamma
        if beta is None:
            beta = default_beta

        self.fft_length = 1 << radix2_exp

        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.bin_per_octave = bin_per_octave
        self.wavelet_type = wavelet_type
        self.scale_type = scale_type
        self.gamma = gamma
        self.beta = beta
        self.thresh = thresh
        self.is_padding = is_padding

        self.order = 1

        fn = self._lib['wsstObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueWSST)),
                       c_int, c_int, POINTER(c_int), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int),
                       POINTER(c_float), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.radix2_exp),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.bin_per_octave)),
           pointer(c_int(self.wavelet_type.value)),
           pointer(c_int(self.scale_type.value)),
           pointer(c_float(self.gamma)),
           pointer(c_float(self.beta)),
           pointer(c_float(self.thresh)),
           pointer(c_int(int(self.is_padding))))

        self._is_created = True

    def get_fre_band_arr(self):
        """
        Get an array of frequency bands of different scales.
        Based on the `scale_type` determination of the initialization.

        Returns
        -------
        out: np.ndarray [shape=(num,)]

        """
        fn = self._lib['wsstObj_getFreBandArr']
        fn.argtypes = [POINTER(OpaqueWSST)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_float * self.num).from_address(p), np.float32).copy()
        return ret

    def get_bin_band_arr(self):
        """
        Get bin band array

        Returns
        -------
        out: np.ndarray [shape=(num,)]
        """
        fn = self._lib['wsstObj_getBinBandArr']
        fn.argtypes = [POINTER(OpaqueWSST)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_int * self.num).from_address(p), np.int32).copy()
        return ret

    def set_order(self, order):
        """
        Set order

        Parameters
        ----------
        order: int, >= 1

        Returns
        -------

        """
        c_fn = self._lib['wsstObj_setOrder']
        c_fn.argtypes = [POINTER(OpaqueWSST), c_int]
        c_fn(self._obj, c_int(order))
        self.order = order

    def wsst(self, data_arr):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., 2**radix2_exp)]
            Input audio data

        Returns
        -------
        m_arr1: np.ndarray [shape=(..., fre, time), dtype=np.complex]
            The matrix of WSST

        m_arr2: np.ndarray [shape=(..., fre, time), dtype=np.complex]
            The matrix of origin(CWT)
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)
        data_arr = check_audio_length(data_arr, self.radix2_exp)

        fn = self._lib['wsstObj_wsst']
        fn.argtypes = [POINTER(OpaqueWSST),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        if data_arr.ndim == 1:
            shape = (self.num, self.fft_length)
            m_real_arr1 = np.zeros(shape, dtype=np.float32)
            m_imag_arr1 = np.zeros(shape, dtype=np.float32)
            m_real_arr2 = np.zeros(shape, dtype=np.float32)
            m_imag_arr2 = np.zeros(shape, dtype=np.float32)

            fn(self._obj, data_arr, m_real_arr1, m_imag_arr1, m_real_arr2, m_imag_arr2)
            m_arr1 = m_real_arr1 + m_imag_arr1 * 1j
            m_arr2 = m_real_arr2 + m_imag_arr2 * 1j
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            shape = (channel_num, self.num, self.fft_length)
            m_real_arr1 = np.zeros(shape, dtype=np.float32)
            m_imag_arr1 = np.zeros(shape, dtype=np.float32)
            m_real_arr2 = np.zeros(shape, dtype=np.float32)
            m_imag_arr2 = np.zeros(shape, dtype=np.float32)

            for i in range(channel_num):
                fn(self._obj, data_arr[i], m_real_arr1[i], m_imag_arr1[i], m_real_arr2[i], m_imag_arr2[i])
            m_arr1 = m_real_arr1 + m_imag_arr1 * 1j
            m_arr2 = m_real_arr2 + m_imag_arr2 * 1j

            m_arr1 = revoke_channel(m_arr1, o_channel_shape, 2)
            m_arr2 = revoke_channel(m_arr2, o_channel_shape, 2)

        m_arr2 = np.ascontiguousarray(m_arr2[..., ::-1, :])
        return m_arr1, m_arr2

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray
        """
        y_coords = self.get_fre_band_arr()
        y_coords = np.insert(y_coords, 0, self.low_fre)
        return y_coords

    def x_coords(self):
        """
        Get the X-axis coordinate

        Returns
        -------
        out: np.ndarray
        """
        x_coords = np.linspace(0, self.fft_length / self.samplate,
                               self.fft_length + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['wsstObj_free']
            free_fn.argtypes = [POINTER(OpaqueWSST)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
