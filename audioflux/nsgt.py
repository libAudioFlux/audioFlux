import warnings
import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.type import (NSGTFilterBankType, SpectralFilterBankScaleType, SpectralFilterBankStyleType,
                            SpectralFilterBankNormalType)
from audioflux.base import Base
from audioflux.utils import check_audio, check_audio_length, note_to_hz

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
    >>> # For radix2_exp=15, then fft_length=2**15=
    >>> audio_arr = audio_arr[:32768]
    array([-5.5879354e-09, -9.3132257e-09,  0.0000000e+00, ...,
           -1.3137090e-01, -1.5649168e-01, -1.8550715e-01], dtype=float32)

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
    array([[1.6688864e+01, 1.6688864e+01, 1.6688864e+01, ..., 3.6051276e+00,
            3.6051276e+00, 3.6051276e+00],
           [1.2269066e+01, 1.2269066e+01, 1.2269066e+01, ..., 2.2744296e+00,
            2.2744296e+00, 2.2744296e+00],
           [8.4105768e+00, 8.4105768e+00, 8.4105768e+00, ..., 3.0433545e+00,
            3.0433545e+00, 3.0433545e+00],
           ...,
           [1.6345095e-02, 1.4781086e-01, 2.2834454e-01, ..., 1.0761657e-02,
            7.4106222e-03, 7.7240127e-03],
           [1.3756215e-02, 5.9653565e-02, 1.4187603e-01, ..., 1.1111894e-03,
            5.4555084e-03, 3.7545679e-04],
           [6.0960581e-03, 4.1618239e-02, 6.6681899e-02, ..., 1.9094696e-03,
            3.1084872e-03, 2.7320804e-03]], dtype=float32)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[0]
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

        self.fft_length = 1 << radix2_exp

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

    def nsgt(self, data_arr):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]
            Input audio data

        Returns
        -------
        m_data_arr: np.ndarray [shape=(fre, time), dtype=np.complex]
            The matrix of NSGT
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)
        data_arr = check_audio_length(data_arr, self.radix2_exp)

        fn = self._lib['nsgtObj_nsgt']
        fn.argtypes = [
            POINTER(OpaqueNSGT),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]

        max_time_length = self.get_max_time_length()
        m_real_arr = np.zeros((self.num, max_time_length), dtype=np.float32)
        m_imag_arr = np.zeros((self.num, max_time_length), dtype=np.float32)

        fn(self._obj, data_arr, m_real_arr, m_imag_arr)
        m_data_arr = m_real_arr + m_imag_arr * 1j

        return m_data_arr

    def get_cell_data(self):
        fn = self._lib['nsgtObj_getCellData']
        fn.argtypes = [
            POINTER(OpaqueNSGT),
            POINTER(POINTER(c_float)),
            POINTER(POINTER(c_float)),
        ]

        pp_m_real_arr = pointer(pointer(c_float()))
        pp_m_imag_arr = pointer(pointer(c_float()))
        fn(self._obj, pp_m_real_arr, pp_m_imag_arr)

        max_time_length = self.get_max_time_length()
        m_real_arr = np.array([pp_m_real_arr.contents[x] for x in range(self.num * max_time_length)], dtype=np.float32)
        m_imag_arr = np.array([pp_m_imag_arr.contents[x] for x in range(self.num * max_time_length)], dtype=np.float32)
        m_real_arr = m_real_arr.reshape(self.num, max_time_length)
        m_imag_arr = m_imag_arr.reshape(self.num, max_time_length)
        return m_real_arr + m_imag_arr * 1j

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
