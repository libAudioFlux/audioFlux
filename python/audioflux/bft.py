from collections import defaultdict
import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.type import (SpectralFilterBankScaleType, SpectralFilterBankStyleType, SpectralFilterBankNormalType,
                            WindowType, SpectralDataType)
from audioflux.base import Base
from audioflux.utils import check_audio, ascontiguous_swapaxex, format_channel, revoke_channel, note_to_hz

__all__ = ["BFT"]


class OpaqueBFT(Structure):
    _fields_ = []


class BFT(Base):
    """
    Based Fourier Transform, similar short-time Fourier transform

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

        Usually set to 12, 24 or 36.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    scale_type: SpectralFilterBankScaleType
        Spectral filter bank type. It determines the type of spectrogram.

        See: `type.SpectralFilterBankScaleType`

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        Linear is not provided.

        See: `type.SpectralFilterBankNormalType`

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    is_reassign: bool
        Whether to use reassign.

    is_temporal: bool
        Whether to get temporal data.

        If True, you can call the **get_temporal_data** method to
        get `energy`/`rms`/`zeroCrossRate` feature.

    See Also
    --------
    NSGT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Create BFT object of Linser(STFT)

    >>> from audioflux.type import (SpectralFilterBankScaleType, SpectralFilterBankStyleType,
    >>>                             WindowType, SpectralDataType)
    >>> obj = af.BFT(num=2049, radix2_exp=12, samplate=sr, low_fre=0., high_fre=16000.,
    >>>              window_type=WindowType.HANN, slide_length=1024,
    >>>              scale_type=SpectralFilterBankScaleType.LINEAR,
    >>>              style_type=SpectralFilterBankStyleType.SLANEY,
    >>>              data_type=SpectralDataType.POWER)

    Extract spectrogram of dB

    >>> import numpy as np
    >>> from audioflux.utils import power_to_db
    >>> spec_arr = obj.bft(audio_arr)
    >>> spec_arr = np.abs(spec_arr)
    >>> spec_dB_arr = power_to_db(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[-1]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='BFT-Linear Spectrogram')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")

    """

    def __init__(self, num, radix2_exp=12, samplate=32000,
                 low_fre=None, high_fre=None, bin_per_octave=12,
                 window_type=WindowType.HANN,
                 slide_length=None,
                 scale_type=SpectralFilterBankScaleType.LINEAR,
                 style_type=SpectralFilterBankStyleType.SLANEY,
                 normal_type=SpectralFilterBankNormalType.NONE,
                 data_type=SpectralDataType.MAG,
                 is_reassign=False, is_temporal=False):
        super(BFT, self).__init__(pointer(OpaqueBFT()))

        self.fft_length = fft_length = 1 << radix2_exp

        # check num
        if num > (fft_length // 2 + 1):
            raise ValueError(f'num={num} is too large')

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

        self.fft_length = fft_length = 1 << radix2_exp
        if slide_length is None:
            slide_length = fft_length // 4

        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.bin_per_octave = bin_per_octave
        self.window_type = window_type
        self.slide_length = slide_length
        self.scale_type = scale_type
        self.style_type = style_type
        self.normal_type = normal_type
        self.data_type = data_type
        self.is_reassign = is_reassign
        self.is_temporal = is_temporal

        self.result_type = 0
        self._temporal_cache = {}
        self._is_temporal_cached = False

        fn = self._lib['bftObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueBFT)), c_int, c_int,
                       POINTER(c_int), POINTER(c_float), POINTER(c_float), POINTER(c_int),
                       POINTER(c_int), POINTER(c_int),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.radix2_exp),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.bin_per_octave)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.scale_type.value)),
           pointer(c_int(self.style_type.value)),
           pointer(c_int(self.normal_type.value)),
           pointer(c_int(self.data_type.value)),
           pointer(c_int(int(self.is_reassign))),
           pointer(c_int(int(self.is_temporal))))
        self._is_created = True

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
        fn = self._lib['bftObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueBFT), c_int]
        return fn(self._obj, c_int(data_length))

    def get_fre_band_arr(self):
        """
        Get an array of frequency bands of different scales.
        Based on the `scale_type` determination of the initialization.

        Returns
        -------
        out: np.ndarray [shape=(fre, )]
        """
        fn = self._lib['bftObj_getFreBandArr']
        fn.argtypes = [POINTER(OpaqueBFT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_float * self.num).from_address(p), np.float32).copy()
        return ret

    def get_bin_band_arr(self):
        """
        Get bin band array

        Returns
        -------
        out: np.ndarray [shape=[n_bin,]]
        """
        fn = self._lib['bftObj_getBinBandArr']
        fn.argtypes = [POINTER(OpaqueBFT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_int * self.num).from_address(p), np.int32).copy()
        return ret

    def set_result_type(self, result_type):
        """
        Set result type.

        Parameters
        ----------
        result_type: int, 0 or 1
            - If `0`, then the result is a matrix of complex numbers.
            - If `1`, then the result is a matrix of real numbers.
        """
        if result_type not in (0, 1):
            raise ValueError(f'`result_type` param error')

        c_fn = self._lib['bftObj_setResultType']
        c_fn.argtypes = [POINTER(OpaqueBFT), c_int]
        c_fn(self._obj, c_int(result_type))
        self.result_type = result_type

    def set_data_norm_value(self, norm_value):
        """
        Set data norm value

        Parameters
        ----------
        norm_value: float
        """
        fn = self._lib['bftObj_setDataNormValue']
        fn.argtypes = [POINTER(OpaqueBFT), c_float]
        fn(self._obj, c_float(norm_value))

    def bft(self, data_arr, result_type=0):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data

        result_type: int，0 or 1
            - If `0`, then the result is a matrix of complex numbers.
            - If `1`, then the result is a matrix of real numbers.

        Returns
        -------
        m_data_arr: np.ndarray [shape=(..., fre, time), dtype=(np.complex or np.float32)]
            The matrix of BFT
        """

        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)
        data_len = data_arr.shape[-1]
        if data_len < self.fft_length:
            raise ValueError(f'radix2_exp={self.radix2_exp}(fft_length={self.fft_length}) '
                             f'is too large for data_arr length={data_len}')
        self._temporal_cache.clear()
        self._is_temporal_cached = False

        if result_type != self.result_type:
            self.set_result_type(result_type)

        fn = self._lib['bftObj_bft']
        fn.argtypes = [POINTER(OpaqueBFT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        time_len = self.cal_time_length(data_len)
        energy_arr, rms_arr, zcr_arr = None, None, None
        if data_arr.ndim == 1:
            m_real_arr = np.zeros((time_len, self.num), dtype=np.float32)
            m_imag_arr = np.zeros((time_len, self.num), dtype=np.float32)
            fn(self._obj, data_arr, c_int(data_len), m_real_arr, m_imag_arr)
            if self.is_temporal:
                energy_arr, rms_arr, zcr_arr = self._get_temporal_data(data_len)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            m_real_arr = np.zeros((channel_num, time_len, self.num), dtype=np.float32)
            m_imag_arr = np.zeros((channel_num, time_len, self.num), dtype=np.float32)
            _temporal_data_dic = defaultdict(list)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], c_int(data_len), m_real_arr[i], m_imag_arr[i])
                if self.is_temporal:
                    _energy_arr, _rms_arr, _zcr_arr = self._get_temporal_data(data_len)
                    _temporal_data_dic['energy'].append(_energy_arr)
                    _temporal_data_dic['rms'].append(_rms_arr)
                    _temporal_data_dic['zcr'].append(_zcr_arr)

            m_real_arr = revoke_channel(m_real_arr, o_channel_shape, 2)
            m_imag_arr = revoke_channel(m_imag_arr, o_channel_shape, 2)
            if self.is_temporal:
                energy_arr = np.stack(_temporal_data_dic['energy'], axis=0)
                rms_arr = np.stack(_temporal_data_dic['rms'], axis=0)
                zcr_arr = np.stack(_temporal_data_dic['zcr'], axis=0)
                energy_arr = revoke_channel(energy_arr, o_channel_shape, 1)
                rms_arr = revoke_channel(rms_arr, o_channel_shape, 1)
                zcr_arr = revoke_channel(zcr_arr, o_channel_shape, 1)

        if self.is_temporal:
            self._temporal_cache['energy'] = energy_arr
            self._temporal_cache['rms'] = rms_arr
            self._temporal_cache['zcr'] = zcr_arr
            self._is_temporal_cached = True
        m_data_arr = (m_real_arr + m_imag_arr * 1j) if self.result_type == 0 else m_real_arr
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        return m_data_arr

    def get_temporal_data(self):
        """
        Get energy/rms/zeroCrossRate feature.

        Need to call `bft` method first.

        Returns
        -------
        energy_arr: np.ndarray [shape=(..., time)]
            energy feature

        rms_arr: np.ndarray [shape=(..., time)]
            rms feature

        zero_cross_arr: np.ndarray [shape=(..., time)]
            zero cross rate feature
        """

        if not self.is_temporal:
            raise ValueError(f'Please set the parameter is_temporal=True when creating the BFT object')
        if not self._is_temporal_cached:
            raise ValueError(f'Please call the `BFT.bft()` method before calling this method')

        energy_arr = self._temporal_cache.get('energy')
        rms_arr = self._temporal_cache.get('rms')
        zcr_arr = self._temporal_cache.get('zcr')
        return energy_arr, rms_arr, zcr_arr

    def _get_temporal_data(self, data_length):
        """
        Call the C function to get energy/rms/zeroCrossRate feature.

        Need to call `bft` method first.

        Parameters
        ----------
        data_length: int
            Audio data length

        Returns
        -------
        energy_arr: np.ndarray [shape=(time,)]
            energy feature

        rms_arr: np.ndarray [shape=(time,)]
            rms feature

        zero_cross_arr: np.ndarray [shape=(time,)]
            zero cross rate feature
        """

        if not self.is_temporal:
            raise ValueError(f'Please set the parameter is_temporal=True when creating the BFT object')

        fn = self._lib['bftObj_getTemporalData']
        fn.argtypes = [POINTER(OpaqueBFT),
                       POINTER(POINTER(c_float)),
                       POINTER(POINTER(c_float)),
                       POINTER(POINTER(c_float)),
                       ]
        time_length = self.cal_time_length(data_length)

        pp_energy_arr = pointer(pointer(c_float()))
        pp_rms_arr = pointer(pointer(c_float()))
        pp_zcr_arr = pointer(pointer(c_float()))

        fn(self._obj,
           pp_energy_arr,
           pp_rms_arr,
           pp_zcr_arr
           )

        energy_arr = np.array([pp_energy_arr.contents[x] for x in range(time_length)], dtype=np.float32)
        rms_arr = np.array([pp_rms_arr.contents[x] for x in range(time_length)], dtype=np.float32)
        zcr_arr = np.array([pp_zcr_arr.contents[x] for x in range(time_length)], dtype=np.float32)

        return energy_arr, rms_arr, zcr_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate

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
        if data_length < self.fft_length:
            raise ValueError(f'radix2_exp={self.radix2_exp}(fft_length={self.fft_length}) '
                             f'is too large for data_length={data_length}')
        x_coords = np.linspace(0, data_length / self.samplate,
                               self.cal_time_length(data_length) + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            fn = self._lib['bftObj_free']
            fn.argtypes = [POINTER(OpaqueBFT)]
            fn.restype = c_void_p
            fn(self._obj)

            self._temporal_cache.clear()
            self._is_temporal_cached = False
