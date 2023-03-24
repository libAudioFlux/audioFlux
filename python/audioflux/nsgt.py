import warnings
import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.type import (NSGTFilterBankType, SpectralFilterBankScaleType, SpectralFilterBankStyleType,
                            SpectralFilterBankNormalType)
from audioflux.base import Base
from audioflux.utils import check_audio, check_audio_length, format_channel, revoke_channel, note_to_hz

__all__ = ['NSGT']


class OpaqueNSGT(Structure):
    _fields_ = []


class NSGT(Base):
    """
    Non-Stationary Gabor Transform (NSGT)

    Parameters
    ----------
    num: int
        Number of frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float or None
        Lowest frequency.

        - Linear/Linsapce/Mel/Bark/Erb, low_fre>=0. `default: 0.0`
        - Octave/Log, low_fre>=32.703. `default: 32.703(C1)`

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

        - Linear is not provided, it is based on ``samplate / (2 ** radix2_exp)``.
        - Octave is not provided, it is based on musical pitch.

    bin_per_octave: int
        Number of bins per octave.

        Only Octave must be provided.

    min_len: int
        Min len

    nsgt_filter_bank_type: NSGTFilterBankType
        NSGT filter bank type.

    scale_type: SpectralFilterBankScaleType
        Spectral filter bank type. It determines the type of spectrogram.

        See: `type.SpectralFilterBankScaleType`

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        The `GAMMATONE` is not supported.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        - Must be set to `NONE` or `BAND_WIDTH`, the `AREA` is not supported.
        - Linear is not provided.

        See: `type.SpectralFilterBankNormalType`

    See Also
    --------
    BFT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)
    >>> # NSGT can only input fft_length data
    >>> # For radix2_exp=15, then fft_length=2**15=32768
    >>> audio_arr = audio_arr[..., :32768]

    Create NSGT object of Octave

    >>> from audioflux.type import (SpectralFilterBankScaleType, SpectralFilterBankStyleType,
    >>>                             SpectralFilterBankNormalType)
    >>> from audioflux.utils import note_to_hz
    >>> obj = af.NSGT(num=84, radix2_exp=15, samplate=sr,
    >>>               low_fre=note_to_hz('C1'), bin_per_octave=12,
    >>>               scale_type=SpectralFilterBankScaleType.OCTAVE,
    >>>               style_type=SpectralFilterBankStyleType.SLANEY,
    >>>               normal_type=SpectralFilterBankNormalType.NONE)

    Extract spectrogram

    >>> import numpy as np
    >>> spec_arr = obj.nsgt(audio_arr)
    >>> spec_arr = np.abs(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[-1]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='NSGT-Octave Spectrogram')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, num=84, radix2_exp=12, samplate=32000,
                 low_fre=None, high_fre=None,
                 bin_per_octave=12, min_len=3,
                 nsgt_filter_bank_type=NSGTFilterBankType.EFFICIENT,
                 scale_type=SpectralFilterBankScaleType.OCTAVE,
                 style_type=SpectralFilterBankStyleType.SLANEY,
                 normal_type=SpectralFilterBankNormalType.BAND_WIDTH):
        super(NSGT, self).__init__(pointer(OpaqueNSGT()))

        self.fft_length = fft_length = 1 << radix2_exp

        # check num
        if num > (fft_length // 2 + 1):
            raise ValueError(f'num={num} is too large')

        # check BPO
        if scale_type == SpectralFilterBankScaleType.OCTAVE and bin_per_octave < 1:
            raise ValueError(f'bin_per_octave={bin_per_octave} must be a positive integer')

        if style_type == SpectralFilterBankStyleType.GAMMATONE:
            raise ValueError(f'style_type={style_type.name} is unsupported')
        if normal_type not in (SpectralFilterBankNormalType.NONE,
                               SpectralFilterBankNormalType.BAND_WIDTH):
            raise ValueError(f'normal_type={normal_type.name} is unsupported')

        if low_fre is None:
            if scale_type in (SpectralFilterBankScaleType.OCTAVE,
                              SpectralFilterBankScaleType.LOG):
                low_fre = note_to_hz('C1')  # 32.703
            else:
                low_fre = 0

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

        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.bin_per_octave = bin_per_octave
        self.min_len = min_len
        self.nsgt_filter_bank_type = nsgt_filter_bank_type
        self.scale_type = scale_type
        self.style_type = style_type
        self.normal_type = normal_type

        fn = self._lib['nsgtObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueNSGT)), c_int, c_int, POINTER(c_int),
                       POINTER(c_float), POINTER(c_float), POINTER(c_int), POINTER(c_int),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)
                       ]
        fn(self._obj,
           c_int(self.num),
           c_int(self.radix2_exp),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.bin_per_octave)),
           pointer(c_int(self.min_len)),
           pointer(c_int(self.nsgt_filter_bank_type.value)),
           pointer(c_int(self.scale_type.value)),
           pointer(c_int(self.style_type.value)),
           pointer(c_int(self.normal_type.value)),
           )
        self._is_created = True

    def get_max_time_length(self):
        """
        Get max time length

        Returns
        -------
        out: int
        """
        fn = self._lib['nsgtObj_getMaxTimeLength']
        fn.argtypes = [POINTER(OpaqueNSGT)]
        fn.restype = c_int
        return fn(self._obj)

    def get_total_time_length(self):
        """
        Get total time length

        Returns
        -------
        out: int
        """
        fn = self._lib['nsgtObj_getTotalTimeLength']
        fn.argtypes = [POINTER(OpaqueNSGT)]
        fn.restype = c_int
        return fn(self._obj)

    def get_time_length_arr(self):
        """
        Get time length array

        Returns
        -------
        out: np.ndarray [shape=(time,)]
        """
        fn = self._lib['nsgtObj_getTimeLengthArr']
        fn.argtypes = [POINTER(OpaqueNSGT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_int * self.num).from_address(p), np.int32).copy()
        return ret

    def get_fre_band_arr(self):
        """
        Get an array of frequency bands of different scales.
        Based on the `scale_type` determination of the initialization.

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        fn = self._lib['nsgtObj_getFreBandArr']
        fn.argtypes = [POINTER(OpaqueNSGT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_float * self.num).from_address(p), np.float32).copy()
        return ret

    def get_bin_band_arr(self):
        """
        Get bin band array

        Returns
        -------
        out: np.ndarray [shape=(n,)]
        """
        fn = self._lib['nsgtObj_getBinBandArr']
        fn.argtypes = [POINTER(OpaqueNSGT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_int * self.num).from_address(p), np.int32).copy()
        return ret

    def set_min_length(self, min_length=3):
        """
        Set min length

        Parameters
        ----------
        min_length: int
        """
        if min_length < 1:
            raise ValueError(f'min_length={min_length} cannot be less than 1')
        fn = self._lib['nsgtObj_setMinLength']
        fn.argtypes = [POINTER(OpaqueNSGT), c_int]
        fn(self._obj, c_int(min_length))
        self.min_len = min_length

    def nsgt(self, data_arr):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., 2**radix2_exp)]
            Input audio data

        Returns
        -------
        m_data_arr: np.ndarray [shape=(..., fre, time), dtype=np.complex]
            The matrix of NSGT
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)
        data_arr = check_audio_length(data_arr, self.radix2_exp)

        fn = self._lib['nsgtObj_nsgt']
        fn.argtypes = [
            POINTER(OpaqueNSGT),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]

        max_time_length = self.get_max_time_length()

        if data_arr.ndim == 1:
            m_real_arr = np.zeros((self.num, max_time_length), dtype=np.float32)
            m_imag_arr = np.zeros((self.num, max_time_length), dtype=np.float32)
            fn(self._obj, data_arr, m_real_arr, m_imag_arr)
            m_nsgt_arr = m_real_arr + m_imag_arr * 1j
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            m_real_arr = np.zeros((channel_num, self.num, max_time_length), dtype=np.float32)
            m_imag_arr = np.zeros((channel_num, self.num, max_time_length), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], m_real_arr[i], m_imag_arr[i])
            m_nsgt_arr = m_real_arr + m_imag_arr * 1j
            m_nsgt_arr = revoke_channel(m_nsgt_arr, o_channel_shape, 2)

        return m_nsgt_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate.

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        y_coords = self.get_fre_band_arr()
        y_coords = np.insert(y_coords, 0, self.low_fre)
        return y_coords

    def x_coords(self, data_length):
        """
        Get the X-axis coordinate

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: np.ndarray [shape=(time,)]
        """
        x_coords = np.linspace(0, data_length * 1. / self.samplate,
                               self.get_max_time_length() + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['nsgtObj_free']
            free_fn.argtypes = [POINTER(OpaqueNSGT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
