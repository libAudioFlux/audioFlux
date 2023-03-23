import ctypes
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p, c_size_t

import numpy as np
from audioflux.base import Base
from audioflux.type import (WindowType, SpectralDataType, SpectralFilterBankType, SpectralFilterBankStyleType,
                            SpectralFilterBankNormalType, ChromaDataNormalType, CepstralRectifyType, CepstralEnergyType,
                            SpectralNoveltyMethodType, SpectralNoveltyDataType)
from audioflux.utils import check_audio, ascontiguous_T, ascontiguous_swapaxex, format_channel, revoke_channel, \
    note_to_hz

__all__ = [
    'Spectrogram',
    'Linear',
    'Mel',
    'Bark',
    'Erb',
    'Chroma',
    'Deep',
    'DeepChroma'
]


class OpaqueSpectrogram(Structure):
    _fields_ = []


class SpectrogramBase(Base):
    def __init__(self, num=0, samplate=32000, low_fre=None, high_fre=None,
                 bin_per_octave=12,
                 radix2_exp=12, window_type=None,
                 slide_length=None,
                 data_type=SpectralDataType.POWER,
                 filter_bank_type=SpectralFilterBankType.LINEAR,
                 filter_style_type=SpectralFilterBankStyleType.SLANEY,
                 filter_normal_type=SpectralFilterBankNormalType.NONE,
                 is_continue=False):
        super(SpectrogramBase, self).__init__(pointer(OpaqueSpectrogram()))

        # check BPO
        if filter_bank_type in (SpectralFilterBankType.OCTAVE, SpectralFilterBankType.OCTAVE_CHROMA):
            if bin_per_octave not in (12, 24, 36):
                raise ValueError(f'bin_per_octave={bin_per_octave} must be 12, 24 or 36')

        if filter_bank_type == SpectralFilterBankType.OCTAVE and num % bin_per_octave != 0:
            raise ValueError(f'num={num} must be an integer multiple of bin_per_octave={bin_per_octave}')

        if low_fre is None:
            if filter_bank_type in (SpectralFilterBankType.LINEAR,
                                    SpectralFilterBankType.LINSPACE,
                                    SpectralFilterBankType.MEL,
                                    SpectralFilterBankType.BARK,
                                    SpectralFilterBankType.ERB):
                low_fre = 0.0
            else:
                low_fre = note_to_hz('C1')  # 32.703

        if high_fre is None:
            high_fre = samplate / 2

        # deep/deep_chroma 默认 window_type=hamm，其他window_type=hann
        if window_type is None:
            if filter_bank_type in (SpectralFilterBankType.DEEP,
                                    SpectralFilterBankType.DEEP_CHROMA):
                window_type = WindowType.HAMM
            else:
                window_type = WindowType.HANN

        # check low_fre
        if filter_bank_type in (SpectralFilterBankType.OCTAVE,
                                SpectralFilterBankType.LOG,
                                SpectralFilterBankType.OCTAVE_CHROMA,
                                SpectralFilterBankType.DEEP,
                                SpectralFilterBankType.DEEP_CHROMA) \
                and low_fre < round(note_to_hz('C1'), 3):
            # octave/log/logChroma/deep/deepChroma>=32.703
            raise ValueError(f'{filter_bank_type.name} low_fre={low_fre} must be greater than or equal to 32.703')
        if low_fre < 0:
            # linear/linspace/mel/bark/erb/chroma low_fre>=0
            raise ValueError(f'{filter_bank_type.name} low_fre={low_fre} must be a non-negative number')

        # check window_type
        if filter_bank_type in (SpectralFilterBankType.DEEP,
                                SpectralFilterBankType.DEEP_CHROMA):
            if window_type not in (WindowType.HAMM, WindowType.HANN, WindowType.RECT):
                raise ValueError(f'{filter_bank_type.name} window_type={window_type.name} must be HAMM, HANN or RECT)')

        fft_length = 1 << radix2_exp
        if slide_length is None:
            slide_length = fft_length // 4

        self.num = num
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.bin_per_octave = bin_per_octave
        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length
        self.is_continue = is_continue
        self.data_type = data_type
        self.filter_bank_type = filter_bank_type
        self.filter_style_type = filter_style_type
        self.filter_normal_type = filter_normal_type

        self.fft_length = fft_length
        self.deep_order = 1

    def set_data_norm_value(self, norm_value):
        """
        Set data norm value

        Parameters
        ----------
        norm_value: float
        """
        fn = self._lib['spectrogramObj_setDataNormValue']
        fn.argtypes = [POINTER(OpaqueSpectrogram), c_float]
        fn(self._obj, c_float(norm_value))

    def set_chroma_data_normal_type(self, data_norm_type):
        """
        Set chroma data normal type

        Parameters
        ----------
        data_norm_type: ChromaDataNormalType
        """
        if not isinstance(data_norm_type, ChromaDataNormalType):
            raise ValueError(f'data_norm_type={data_norm_type} must be of type')

        fn = self._lib['spectrogramObj_setChromaDataNormalType']
        fn.argtypes = [POINTER(OpaqueSpectrogram), c_int]
        fn(self._obj, c_int(data_norm_type.value))

    def set_deep_order(self, deep_order):
        """
        Set deep order

        Parameters
        ----------
        deep_order: int [1,2,3,4]
            spectrogram shape
                | If deep_order=1/2 -> 3*timeLength*num
                | If deep_order=3/4 -> 5*timeLength*num

        """
        if deep_order not in (1, 2, 3, 4):
            raise ValueError(f'deep_order={deep_order} must be in (1, 2, 3, 4)')

        fn = self._lib['spectrogramObj_setDeepOrder']
        fn.argtypes = [POINTER(OpaqueSpectrogram), c_int]
        self.deep_order = deep_order
        fn(self._obj, c_int(deep_order))

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
            time length
        """
        fn = self._lib['spectrogramObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueSpectrogram), c_int]
        return fn(self._obj, c_int(data_length))

    def get_fre_band_arr(self):
        """
        Get frequency band array

        Returns
        -------
        out: np.ndarray [shape=(num,)]

        """
        fn = self._lib['spectrogramObj_getFreBandArr']
        fn.argtypes = [POINTER(OpaqueSpectrogram)]
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
        fn = self._lib['spectrogramObj_getBinBandArr']
        fn.argtypes = [POINTER(OpaqueSpectrogram)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_int * self.num).from_address(p), np.int32).copy()
        return ret

    def get_band_num(self):
        """
        Get band num

        Returns
        -------
        out: int
        """
        fn = self._lib['spectrogramObj_getBandNum']
        fn.argtypes = [POINTER(OpaqueSpectrogram)]
        fn.restype = c_int
        ret = fn(self._obj)
        return ret

    def get_bin_band_length(self):
        """
        Get bin band length

        Returns
        -------
        out: int
        """
        fn = self._lib['spectrogramObj_getBinBandLength']
        fn.argtypes = [POINTER(OpaqueSpectrogram)]
        fn.restype = c_int
        ret = fn(self._obj)
        return ret

    def spectrogram(self, data_arr, is_phase_arr=False):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data
        is_phase_arr: bool
            Whether to return phase array

        Returns
        -------
        out: np.ndarray [shape=(..., fre, time) or (..., deep, fre, time)]
            Spectrogram matrix array.

            If filter_bank_type is ``DEEP`` or ``DEEP_CHROMA``, the shape is ``(..., deep, fre, time)``;
            otherwise the shep is ``(..., fre, time)``.
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)
        data_len = data_arr.shape[-1]
        if data_len < self.fft_length:
            raise ValueError(f'radix2_exp={self.radix2_exp}(fft_length={self.fft_length}) '
                             f'is too large for data_arr length={data_len}')

        if is_phase_arr and self.filter_bank_type is not SpectralFilterBankType.LINEAR:
            raise ValueError(f'Only LINEAR bank type has phase arr')

        # deep type has 3 dimensions, other has 2 dimensions
        spec_arr_dim_num = 3 if self.filter_bank_type == SpectralFilterBankType.DEEP else 2

        fn = self._lib['spectrogramObj_spectrogram']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=spec_arr_dim_num, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=spec_arr_dim_num,
                                   flags='C_CONTIGUOUS') if is_phase_arr else POINTER(c_int)
        ]

        time_len = self.cal_time_length(data_len)

        # deep bank need to set 3 dimensions
        if spec_arr_dim_num == 3:
            # 1/2->3*timeLength*num
            # 3/4->5*timeLength*num
            if self.deep_order in (1, 2):
                deep_len = 3
            else:
                deep_len = 5
            size = (deep_len, time_len, self.num)
        else:
            size = (time_len, self.num)

        if data_arr.ndim == 1:
            m_spec_arr = np.zeros(size, dtype=np.float32)
            m_phase_arr = None
            if is_phase_arr:
                m_phase_arr = np.zeros((time_len, int((1 << self.radix2_exp) / 2 + 1)), dtype=np.float32)

            fn(self._obj, data_arr, c_int(data_len), m_spec_arr, m_phase_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            size = (channel_num, *size)
            m_spec_arr = np.zeros(size, dtype=np.float32)
            m_phase_arr = None
            if is_phase_arr:
                m_phase_arr = np.zeros((channel_num, time_len, int((1 << self.radix2_exp) / 2 + 1)), dtype=np.float32)

            for i in range(channel_num):
                _m_phase_arr = m_phase_arr[i] if is_phase_arr else None
                fn(self._obj, data_arr[i], c_int(data_len), m_spec_arr[i], _m_phase_arr)

            m_spec_arr = revoke_channel(m_spec_arr, o_channel_shape, spec_arr_dim_num)
            if is_phase_arr:
                m_phase_arr = revoke_channel(m_phase_arr, o_channel_shape, 2)

        m_spec_arr = ascontiguous_swapaxex(m_spec_arr, -1, -2)

        if is_phase_arr:
            m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)
            return m_spec_arr, m_phase_arr
        else:
            return m_spec_arr

    def deconv(self, m_data_arr):
        """
        Compute the spectral deconv feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        m_tone_arr: np.ndarray [shape=(..., time)]
            The matrix of tone
        m_pitch_arr: np.ndarray [shape=(..., time)]
            The matrix of pitch
        """
        if m_data_arr.ndim != 2:
            raise ValueError(f'm_data_arr must be ')
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_deconv']
        fn.argtypes = [POINTER(OpaqueSpectrogram),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        shape = m_data_arr.shape  # (time, fre)
        m_tone_arr = np.zeros(shape, dtype=np.float32)
        m_pitch_arr = np.zeros(shape, dtype=np.float32)
        fn(self._obj, m_data_arr, m_tone_arr, m_pitch_arr)

        m_tone_arr = ascontiguous_T(m_tone_arr)
        m_pitch_arr = ascontiguous_T(m_pitch_arr)
        return m_tone_arr, m_pitch_arr

    def mfcc(self, m_data_arr, cc_num=13):
        """
        Mel-frequency cepstral coefficients (MFCC)

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data
        cc_num: int
            mfcc num, usually set to 13, 20 or 40

        Returns
        -------
        out: np.ndarray [shape=(cc_num, time)]
        """
        if not (self.filter_bank_type == SpectralFilterBankType.MEL
                and self.filter_style_type == SpectralFilterBankStyleType.SLANEY):
            raise ValueError(f'``filter_bank_type`` must be ``MEL`` and '
                             f'``filter_style_type`` must be ``SLANEY``')
        if cc_num > self.num:
            raise ValueError(f'cc_num={cc_num} must be less than num={self.num}')
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_mfcc']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]

        time_len, _ = m_data_arr.shape
        ret = np.zeros((time_len, cc_num), dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(cc_num), ret)
        return ascontiguous_T(ret)

    def bfcc(self, m_data_arr, cc_num=13):
        """
        Bark-frequency cepstral coefficients (BFCC)

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data
        cc_num: int
            bfcc num, usually set to 13, 20 or 40

        Returns
        -------
        out: np.ndarray [shape=(cc_num, time)]
        """

        if not (self.filter_bank_type == SpectralFilterBankType.BARK
                and self.filter_style_type == SpectralFilterBankStyleType.SLANEY):
            raise ValueError(f'``filter_bank_type`` must be ``BARK`` and '
                             f'``filter_style_type`` must be ``SLANEY``')
        if cc_num > self.num:
            raise ValueError(f'cc_num={cc_num} must be less than num={self.num}')
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_bfcc']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]

        time_len, _ = m_data_arr.shape
        ret = np.zeros((time_len, cc_num), dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(cc_num), ret)
        return ascontiguous_T(ret)

    def gtcc(self, m_data_arr, cc_num=13):
        """
        Gammatone cepstral coefficients (GTCC)

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data
        cc_num: int
            gtcc num, usually set to 13, 20 or 40

        Returns
        -------
        out: np.ndarray [shape=(cc_num, time)]
        """
        if not (self.filter_bank_type == SpectralFilterBankType.ERB
                and self.filter_style_type == SpectralFilterBankStyleType.GAMMATONE):
            raise ValueError(f'``filter_bank_type`` must be ``ERB`` and '
                             f'``filter_style_type`` must be ``GAMMATONE``')
        if cc_num > self.num:
            raise ValueError(f'cc_num={cc_num} must be less than num={self.num}')
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_gtcc']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]
        time_len, _ = m_data_arr.shape
        ret = np.zeros((time_len, cc_num), dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(cc_num), ret)
        return ascontiguous_T(ret)

    def xxcc(self, m_data_arr, cc_num=13, rectify_type=CepstralRectifyType.LOG):
        """
        xx cepstral coefficients

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data
        cc_num: int
            xxcc num, usually set to 13, 20 or 40
        rectify_type: CepstralRectifyType

        Returns
        -------
        out: np.ndarray [shape=(cc_num, time)]
        """
        if cc_num > self.num:
            raise ValueError(f'cc_num={cc_num} must be less than num={self.num}')
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_xxcc']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')
        ]
        time_len, _ = m_data_arr.shape
        ret = np.zeros((time_len, cc_num), dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(cc_num), pointer(c_int(rectify_type.value)), ret)
        return ascontiguous_T(ret)

    def set_edge(self, start: int, end: int):
        """
        Set edge

        Parameters
        ----------
        start: int
            0 ~ end
        end: int
            start ~ num-1

        Returns
        -------

        """
        if not 0 <= start < end:
            raise ValueError(f'start={start} must be in range [0, {end})')
        if not start < end <= self.num - 1:
            raise ValueError(f'start={end} must be in range ({start}, {self.num - 1}]')

        fn = self._lib['spectrogramObj_setEdge']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            c_int,
            c_int,
        ]
        fn(self._obj, c_int(start), c_int(end))

    def set_edge_arr(self, index_arr):
        """
        Set edge arr

        Parameters
        ----------
        index_arr: np.ndarray [shape=(num,)] or list
            fre下标

        Returns
        -------

        """

        fn = self._lib['spectrogramObj_setEdgeArr']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            POINTER(c_int),
            c_int
        ]

        index_len = len(index_arr)

        calloc_fn = self._lib['calloc']
        calloc_fn.argtypes = [c_size_t, c_size_t]
        calloc_fn.restype = c_void_p
        address = calloc_fn(c_size_t(index_len), c_size_t(ctypes.sizeof(c_int)))

        p_index_arr = (c_int * index_len).from_address(address)
        for i, v in enumerate(index_arr):
            p_index_arr[i] = v
        p_index_arr = ctypes.cast(p_index_arr, POINTER(c_int))

        fn(self._obj, p_index_arr, c_int(index_len))

    def preprocess(self, m_data_arr_1, m_data_arr_3=None):
        """
        Pre-processing spectrum data

        You need to execute before calling the spectral correlation function

        Parameters
        ----------
        m_data_arr_1: np.ndarray [shape=(fre, time)]
            Spectrogram data
        m_data_arr_3: np.ndarray or None
        """
        m_data_arr_1 = ascontiguous_T(m_data_arr_1)

        fn = self._lib['spectrogramObj_preprocess']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            POINTER(c_void_p) if m_data_arr_3 is None else np.ctypeslib.ndpointer(dtype=np.float32,
                                                                                  ndim=2,
                                                                                  flags='C_CONTIGUOUS')
        ]
        fn(self._obj, m_data_arr_1, m_data_arr_3)

    def flatness(self, m_data_arr):
        """
        Compute the spectral flatness feature.

        :math:`\qquad flatness=\\frac{\\left ( \prod_{k=b_1}^{b_2} s_k  \\right)^{ \\frac{1}{b_2-b_1} } } {\\frac{1}{b_2-b_1} \sum_{ k=b_1 }^{b_2} s_k}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        flatness: np.ndarray [shape=(time,)]
            flatness frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_flatness']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * f
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def flux(self, m_data_arr, step=1, p=2, is_positive=False, is_no_exp=True, tp=0):
        """
        Compute the spectral flux feature.

        :math:`\qquad flux(t)=\\left( \sum_{k=b_1}^{b_2} |s_k(t)-s_k(t-1) |^{p}  \\right)^{\\frac{1}{p}}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        * In general :math:`s_k(t) \geq s_k(t-1)` participate in the calculation

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        step: int
            Compute time axis steps, like 1/2/3/...

        p: int, 1 or 2
            norm: 1 abs; 2 pow

        is_positive: bool
            Whether to set negative numbers to 0

        is_no_exp: bool
            Whether to exp

        tp: int, 0 or 1
            0 sum 1 mean

        Returns
        -------
        flux: np.ndarray [shape=(time,)]
            flux frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_flux']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int, c_float, c_int,
            POINTER(c_int), POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        time_len, _ = m_data_arr.shape  # t * n
        ret = np.zeros(time_len, dtype=np.float32)
        fn(self._obj,
           m_data_arr,
           c_int(step),
           c_float(p),
           c_int(int(is_positive)),
           pointer(c_int(int(is_no_exp))),
           pointer(c_int(tp)),
           ret)
        return ret

    def rolloff(self, m_data_arr, threshold=0.95):
        """
        Compute the spectral rolloff feature.

        :math:`\qquad \sum_{k=b_1}^{i}|s_k| \geq \eta \sum_{k=b_1}^{b_2}s_k`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`\eta \in (0,1)`, generally take 0.95 or 0.85, satisfy the condition :math:`i` get :math:`f_i` rolloff frequency

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        threshold: float, [0,1]
            rolloff threshold. Generally take 0.95 or 0.85.

        Returns
        -------
        rolloff: np.ndarray [shape=(time,)]
            rolloff frequency for each time
        """
        if not 0 <= threshold <= 1:
            raise ValueError(f'threshold is between 0 and 1')

        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_rolloff']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_float(threshold), ret)
        return ret

    def centroid(self, m_data_arr):
        """
        Compute the spectral centroid feature.

        :math:`\qquad \mu_1=\\frac{\sum_{ k=b_1 }^{b_2} f_ks_k } {\sum_{k=b_1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        spectral: np.ndarray [shape=(time,)]
            spectral frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_centroid']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def spread(self, m_data_arr):
        """
        Compute the spectral spread feature.

        :math:`\qquad \mu_2=\sqrt{\\frac{\sum_{ k=b_1 }^{b_2} (f_k-\mu_1)^2 s_k } {\sum_{k=b_1}^{b_2} s_k } }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`u_1`: `Spectral.centroid`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        spread: np.ndarray [shape=(time,)]
            spread frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_spread']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def skewness(self, m_data_arr):
        """
        Compute the spectral skewness feature.

        :math:`\qquad \mu_3=\\frac{\sum_{ k=b_1 }^{b_2} (f_k-\mu_1)^3 s_k } {(\mu_2)^3 \sum_{k=b_1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`u_1`: `Spectral.centroid`
        * :math:`u_2`: `Spectral.spread`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        skewness: np.ndarray [shape=(time,)]
            skewness frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_skewness']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def kurtosis(self, m_data_arr):
        """
        Compute the spectral kurtosis feature.

        :math:`\qquad \mu_4=\\frac{\sum_{ k=b_1 }^{b_2} (f_k-\mu_1)^4 s_k } {(\mu_2)^4 \sum_{k=b_1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`u_1`: `Spectral.centroid`
        * :math:`u_2`: `Spectral.spread`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        kurtosis: np.ndarray [shape=(time,)]
            kurtosis frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_kurtosis']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def entropy(self, m_data_arr, is_norm=False):
        """
        Compute the spectral entropy feature.

        Set: :math:`p_k=\\frac{s_k}{\sum_{k=b_1}^{b_2}s_k}`

        :math:`\qquad entropy1= \\frac{-\sum_{ k=b_1 }^{b_2} p_k \log(p_k)} {\log(b_2-b_1)}`

        Or

        :math:`\qquad entropy2= {-\sum_{ k=b_1 }^{b_2} p_k \log(p_k)}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        is_norm: bool
            Whether to norm

        Returns
        -------
        entropy: np.ndarray [shape=(time,)]
            entropy frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_entropy']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(int(is_norm)), ret)
        return ret

    def crest(self, m_data_arr):
        """
        Compute the spectral crest feature.

        :math:`\qquad crest =\\frac{max(s_{k\in_{[b_1,b_2]} }) } {\\frac{1}{b_2-b_1} \sum_{ k=b_1 }^{b_2} s_k}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        crest: np.ndarray [shape=(time,)]
            crest frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_crest']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def slope(self, m_data_arr):
        """
        Compute the spectral slope feature.

        :math:`\qquad slope=\\frac{ \sum_{k=b_1}^{b_2}(f_k-\mu_f)(s_k-\mu_s) } { \sum_{k=b_1}^{b_2}(f_k-\mu_f)^2 }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`\mu_f`: average frequency value
        * :math:`\mu_s`: average spectral value

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        slope: np.ndarray [shape=(time,)]
            slope frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_slope']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def decrease(self, m_data_arr):
        """
        Compute the spectral decrease feature.

        :math:`\qquad decrease=\\frac { \sum_{k=b_1+1}^{b_2} \\frac {s_k-s_{b_1}}{k-1} } { \sum_{k=b_1+1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        decrease: np.ndarray [shape=(time,)]
            decrease frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_decrease']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def band_width(self, m_data_arr, p=2):
        """
        Compute the spectral band_width feature.

        :math:`\qquad bandwidth=\\left(\sum_{k=b_1}^{b_2} s_k(f_k-centroid)^p \\right)^{\\frac{1}{p}}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * centroid: `Spectral.centroid`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        p: int, 1 or 2
            norm: 1 abs; 2 pow

        Returns
        -------
        band_width: np.ndarray [shape=(time,)]
            band_width frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_bandWidth']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_float(p), ret)
        return ret

    def rms(self, m_data_arr):
        """
        Compute the spectral rms feature.

        :math:`\qquad rms=\sqrt{ \\frac{1}{N} \sum_{n=1}^N x^2[n] }=\sqrt {\\frac{1}{N^2}\sum_{m=1}^N |X[m]|^2 }`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        rms: np.ndarray [shape=(time,)]
            rms frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_rms']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def energy(self, m_data_arr, is_log=False, gamma=10.):
        """
        Compute the spectral energy feature.

        :math:`\qquad energy=\sum_{n=1}^N x^2[n] =\\frac{1}{N}\sum_{m=1}^N |X[m]|^2`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        is_log: bool
            Whether to log

        gamma: float
            energy gamma value.

        Returns
        -------
        energy: np.ndarray [shape=(time,)]
            energy frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_energy']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, int(is_log), gamma, ret)
        return ret

    def hfc(self, m_data_arr):
        """
        Compute the spectral hfc feature.

        :math:`\qquad hfc(t)=\\frac{\sum_{k=b_1}^{b_2} s_k(t)k }{b_2-b_1+1}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        hfc: np.ndarray [shape=(time,)]
            hfc frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_hfc']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, ret)
        return ret

    def sd(self, m_data_arr, step=1, is_positive=False):
        """
        Compute the spectral sd feature.

        :math:`\qquad sd(t)=flux(t)`

        satisfies the calculation of :math:`s_k(t) \ge s_k(t-1)`, :math:`p=2`，the result is not :math:`1/p`

        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * flux: `Spectral.flux`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        step: int
            Compute time axis steps, like 1/2/3/...

        is_positive: bool
            Whether to set negative numbers to 0

        Returns
        -------
        sd: np.ndarray [shape=(time,)]
            sd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_sd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(step), c_int(int(is_positive)), ret)
        return ret

    def sf(self, m_data_arr, step=1, is_positive=False):
        """
        Compute the spectral sf feature.

        :math:`\qquad sf(t)=flux(t)`

        satisfies the calculation of :math:`s_k(t) \ge s_k(t-1)`, :math:`p=1`

        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * flux: `Spectral.flux`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        step: int
            Compute time axis steps, like 1/2/3/...

        is_positive: bool
            Whether to set negative numbers to 0

        Returns
        -------
        sf: np.ndarray [shape=(time,)]
            sf frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_sf']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(step), c_int(int(is_positive)), ret)
        return ret

    def mkl(self, m_data_arr, tp=0):
        """
        Compute the spectral mkl feature.

        :math:`\qquad mkl(t)=\sum_{k=b_1}^{b_2} \log\\left(1+ \cfrac {s_k(t)}{s_k(t-1)} \\right)`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        tp: int, 0 or 1
            0 sum 1 mean

        Returns
        -------
        mkl: np.ndarray [shape=(time,)]
            mkl frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_mkl']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(tp), ret)
        return ret

    def pd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral pd feature.

        :math:`\qquad \psi_k(t)` is set as the phase function of point `k` at time `t`.

        :math:`\qquad \psi_k^{\prime}(t)=\psi_k(t)-\psi_k(t-1)`

        :math:`\qquad \psi_k^{\prime\prime}(t)=\psi_k^{\prime}(t)-\psi_k^{\prime}(t-1) = \psi_k(t)-2\psi_k(t-1)+\psi_k(t-2)`

        :math:`\qquad pd(t)= \\frac {\sum_{k=b_1}^{b_2} \|  \psi_k^{\prime\prime}(t) \|} {b_2-b_1+1}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        m_phase_arr: np.ndarray [shape=(fre, time)]
            Phase data

        Returns
        -------
        pd: np.ndarray [shape=(time,)]
            pd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)
        m_phase_arr = ascontiguous_T(m_phase_arr)

        fn = self._lib['spectrogramObj_pd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, m_phase_arr, ret)
        return ret

    def wpd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral wpd feature.

        :math:`\qquad \psi_k(t)` is set as the phase function of point k at time t.

        :math:`\qquad \psi_k^{\prime}(t)=\psi_k(t)-\psi_k(t-1)`

        :math:`\qquad \psi_k^{\prime\prime}(t)=\psi_k^{\prime}(t)-\psi_k^{\prime}(t-1) = \psi_k(t)-2\psi_k(t-1)+\psi_k(t-2)`

        :math:`\qquad wpd(t)= \\frac {\sum_{k=b_1}^{b_2}  \| \psi_k^{\prime\prime}(t) \|s_k(t)}{b_2-b_1+1}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        m_phase_arr: np.ndarray [shape=(fre, time)]
            Phase data

        Returns
        -------
        wpd: np.ndarray [shape=(time,)]
            wpd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)
        m_phase_arr = ascontiguous_T(m_phase_arr)

        fn = self._lib['spectrogramObj_wpd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, m_phase_arr, ret)
        return ret

    def nwpd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral nwpd feature.

        :math:`\qquad nwpd(t)= \\frac {wpd} {\mu_s}`

        * wpd: `Spectral.wpd`
        * :math:`\mu_s`: the mean of :math:`s_k(t)`
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        m_phase_arr: np.ndarray [shape=(fre, time)]
            Phase data

        Returns
        -------
        nwpd: np.ndarray [shape=(time,)]
            nwpd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)
        m_phase_arr = ascontiguous_T(m_phase_arr)

        fn = self._lib['spectrogramObj_nwpd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, m_phase_arr, ret)
        return ret

    def cd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral cd feature.

        :math:`\qquad \psi_k(t)` is set as the phase function of point k at time t.

        :math:`\qquad \\alpha_k(t)=s_k(t) e^{j(2\psi_k(t)-\psi_k(t-1))}`

        :math:`\qquad \\beta_k(t)=s_k(t) e^{j\psi_k(t)}`

        :math:`\qquad cd(t)=\sum_{k=b_1}^{b_2} \| \\beta_k(t)-\\alpha_k(t-1)  \|`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        m_phase_arr: np.ndarray [shape=(fre, time)]
            Phase data

        Returns
        -------
        cd: np.ndarray [shape=(time,)]
            cd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)
        m_phase_arr = ascontiguous_T(m_phase_arr)

        fn = self._lib['spectrogramObj_cd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, m_phase_arr, ret)
        return ret

    def rcd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral rcd feature.

        :math:`\qquad rcd(t)=cd`

        participate in the sum calculation when :math:`s_k(t) \geq s_k(t-1)` is satisfied

        * cd: `Spectral.cd`
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        m_phase_arr: np.ndarray [shape=(fre, time)]
            Phase data

        Returns
        -------
        rcd: np.ndarray [shape=(time,)]
            rcd frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)
        m_phase_arr = ascontiguous_T(m_phase_arr)

        fn = self._lib['spectrogramObj_rcd']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, m_phase_arr, ret)
        return ret

    def broadband(self, m_data_arr, threshold):
        """
        Compute the spectral broadband feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        threshold: float, [0,1]
            broadband threshold

        Returns
        -------
        broadband: np.ndarray [shape=(time,)]
            broadband frequency for each time
        """
        if not 0 <= threshold <= 1:
            raise ValueError(f'threshold must be 0 or 1')

        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_broadband']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_float(threshold), ret_arr)
        return ret_arr

    def novelty(self, m_data_arr, step=1, threshold=0.,
                method_type=SpectralNoveltyMethodType.SUB,
                data_type=SpectralNoveltyDataType.VALUE):
        """
        Compute the spectral novelty feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        step: int
            Compute time axis steps, like 1/2/3/...

        threshold: float [0,1]
            Novelty threshold.

        method_type: SpectralNoveltyMethodType
            Novelty method type.

        data_type: SpectralNoveltyDataType
            Novelty data type.

        Returns
        -------
        novelty: np.ndarray [shape=(time,)]
            Novelty frequency per time step.
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_novelty']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            POINTER(c_int),
            POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(step), c_float(threshold),
           pointer(c_int(method_type.value)),
           pointer(c_int(data_type.value)),
           ret_arr)
        return ret_arr

    def eef(self, m_data_arr, is_norm=False):
        """
        Compute the spectral eef feature.

        :math:`\qquad p_k=\\frac{s_k}{\sum_{k=b_1}^{b_2}s_k}`

        :math:`\qquad entropy2= {-\sum_{ k=b_1 }^{b_2} p_k \log(p_k)}`

        :math:`\qquad eef=\sqrt{ 1+| energy\\times entropy2| }`

        * energy: `Spectral.energy`
        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        is_norm: bool
            Whether to norm

        Returns
        -------
        eef: np.ndarray [shape=(time,)]
            eef frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_eef']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(int(is_norm)), ret_arr)
        return ret_arr

    def eer(self, m_data_arr, is_norm=False, gamma=1.):
        """
        Compute the spectral eer feature.

        :math:`\qquad le=\log_{10}(1+\gamma \\times energy), \gamma \in (0,\infty)`, represents `log` compression of data

        :math:`\qquad p_k=\\frac{s_k}{\sum_{k=b_1}^{b_2}s_k}`

        :math:`\qquad entropy2= {-\sum_{ k=b_1 }^{b_2} p_k \log(p_k)}`

        :math:`\qquad eer=\\sqrt{ 1+\\left| \\cfrac{le}{entropy2}\\right| }`

        * energy: Spectral.energy
        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        is_norm: bool
            Whether to norm

        gamma: float
            Usually set is 1./10./20.etc, song is 0.5

        Returns
        -------
        eer: np.ndarray [shape=(time,)]
            eer frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_eer']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        ret_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, c_int(int(is_norm)), c_float(gamma), ret_arr)
        return ret_arr

    def max(self, m_data_arr):
        """
        Compute the spectral max feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        val_arr: np.ndarray [shape=(time,)]
            max value for each time

        fre_arr: np.ndarray [shape=(time,)]
            max frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_max']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        val_arr = np.zeros(n_len, dtype=np.float32)
        fre_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, val_arr, fre_arr)
        return val_arr, fre_arr

    def mean(self, m_data_arr):
        """
        Compute the spectral mean feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        val_arr: np.ndarray [shape=(time,)]
            mean value for each time

        fre_arr: np.ndarray [shape=(time,)]
            mean frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_mean']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        val_arr = np.zeros(n_len, dtype=np.float32)
        fre_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, val_arr, fre_arr)
        return val_arr, fre_arr

    def var(self, m_data_arr):
        """
        Compute the spectral var feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            Spectrogram data

        Returns
        -------
        val_arr: np.ndarray [shape=(time,)]
            var value for each time

        fre_arr: np.ndarray [shape=(time,)]
            var frequency for each time
        """
        m_data_arr = ascontiguous_T(m_data_arr)

        fn = self._lib['spectrogramObj_var']
        fn.argtypes = [
            POINTER(OpaqueSpectrogram),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        n_len, m_len = m_data_arr.shape  # t * n
        val_arr = np.zeros(n_len, dtype=np.float32)
        fre_arr = np.zeros(n_len, dtype=np.float32)
        fn(self._obj, m_data_arr, val_arr, fre_arr)
        return val_arr, fre_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        if self.filter_bank_type in (SpectralFilterBankType.CHROMA,
                                     SpectralFilterBankType.OCTAVE_CHROMA,
                                     SpectralFilterBankType.DEEP_CHROMA):
            y_coords = np.arange(0, self.num + 1)
        else:
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
        x_coords = np.arange(0, self.cal_time_length(data_length) + 1) * (self.slide_length / self.samplate)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['spectrogramObj_free']
            free_fn.argtypes = [POINTER(OpaqueSpectrogram)]
            free_fn.restype = c_void_p
            free_fn(self._obj)


class Spectrogram(SpectrogramBase):
    """
    Generic spectrogram class.

    It can create spectrogram objects of
    Linear/Linspace/Mel/Bark/Erb/Octave/Log/Deep/Chroma/LogChroma/DeepChroma

    Parameters
    ----------
    num: int
        Number of frequency bins to generate, starting at `low_fre`.

        - Linear: not provided. default: ``num = fftLength / 2 + 1``
        - Linspace: Similar to Linear, but you can set it at will
        - Mel/Bark/Erb: Usuall use 128/256/...
        - Octave: Based on pitch. ``num=octave*bin_per_octave``
        - Log: Similar to Octave, but you can set it at will
        - Deep: Based on pitch, but you can set it at will
        - Chroma/LogChroma/DeepChroma: Based on pitch, usuall use 12

    samplate: int
        Sampling rate of the incoming audio

    low_fre: float or None
        Lowest frequency

        - Linear/Linsapce/Mel/Bark/Erb/Chroma, low_fre>=0. `default: 0.0`
        - Octave/Log/LogChroma/Deep/DeepChroma, low_fre>=32.703. `default: 32.703(C1)`

    high_fre: float
        Highest frequency

        Only Linear/Linspace/Mel/Bark/Erb/Log/Chroma/LogChroma must be provided. `default: 16000(samplate/2)`

    bin_per_octave: int
        Number of bins per octave.

        Only Octave/LogChroma must be provided.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    window_type: WindowType
        Window type

        - The default value of `Deep/DeepChroma` is **hamm**, and others are **hann**
        - Deep/DeepChroma can only be set to `hamm/hann/rect`

        See: `type.WindowType`

    slide_length: int
        Window sliding length

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type, you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    filter_bank_type: SpectralFilterBankType
        Spectral filter bank type. It determines the type of spectrogram, like Linear/Mel/Bank/Erb etc.

    filter_style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    filter_normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        Linear/Chroma/Deep/DeepChroma is not provided.

        See: `type.SpectralFilterBankNormalType`

    Notes
    -----
    * Chroma
        * No spectral related features
        * Does not exist freBandArr/binBandArr

    See Also
    --------
    Mel
    Bark
    Erb
    Chroma
    Deep
    DeepChroma

    Examples
    --------
    .. plot::

        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get Linear Spectrogram data

        >>> from audioflux.spectrogram import Spectrogram
        >>> from audioflux.type import SpectralFilterBankType
        >>> spec_obj = Spectrogram(samplate=32000, radix2_exp=12,
        ...                        filter_bank_type=SpectralFilterBankType.LINEAR)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  y_coords=spec_obj.y_coords(),
        ...                  scale='log', title='Linear')


    """

    def __init__(self, num=0, samplate=32000, low_fre=None, high_fre=None,
                 bin_per_octave=12,
                 radix2_exp=12, window_type=None,
                 slide_length=None,
                 data_type=SpectralDataType.POWER,
                 filter_bank_type=SpectralFilterBankType.LINEAR,
                 filter_style_type=SpectralFilterBankStyleType.SLANEY,
                 filter_normal_type=SpectralFilterBankNormalType.NONE):
        super(Spectrogram, self).__init__(num=num, samplate=samplate,
                                          low_fre=low_fre, high_fre=high_fre,
                                          bin_per_octave=bin_per_octave,
                                          radix2_exp=radix2_exp, window_type=window_type,
                                          slide_length=slide_length,
                                          data_type=data_type,
                                          filter_bank_type=filter_bank_type,
                                          filter_style_type=filter_style_type,
                                          filter_normal_type=filter_normal_type)

        fn = self._lib['spectrogramObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)),  # obj
                       c_int,  # num
                       POINTER(c_int),  # samplate
                       POINTER(c_float),  # low_fre
                       POINTER(c_float),  # high_fre
                       POINTER(c_int),  # bin_per_octave
                       POINTER(c_int),  # radix2_exp
                       POINTER(c_int),  # window_type
                       POINTER(c_int),  # slide_length
                       POINTER(c_int),  # is_continue
                       POINTER(c_int),  # data_type
                       POINTER(c_int),  # filter_bank_type
                       POINTER(c_int),  # filter_style_type
                       POINTER(c_int)]  # filter_normal_type
        fn(self._obj,
           c_int(self.num),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.bin_per_octave)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(int(self.is_continue))),
           pointer(c_int(self.data_type.value)),
           pointer(c_int(self.filter_bank_type.value)),
           pointer(c_int(self.filter_style_type.value)),
           pointer(c_int(self.filter_normal_type.value)))

        # 计算linear的num
        if filter_bank_type == SpectralFilterBankType.LINEAR:
            self.num = self.get_band_num()

        self._is_created = True


class Linear(SpectrogramBase):
    """
    Simple linear spectrogram class.

    It can create simple linear spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Linear class fixed parameter:
        * fre range: 0 - 16000(samplate/2)
        * slide_length: ``fft_length // 4``
        * window_type: HANN
        * data_type: POWER
        * style_type: SLANEY

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::

        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get Linear Spectrogram data

        >>> from audioflux.spectrogram import Linear
        >>> spec_obj = Linear(samplate=32000, radix2_exp=12)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  y_coords=spec_obj.y_coords(),
        ...                  scale='log', title='Linear')
    """

    def __init__(self, samplate=32000, radix2_exp=12):
        super(Linear, self).__init__(num=0, samplate=samplate, low_fre=0.0, high_fre=16000.0,
                                     bin_per_octave=12,
                                     radix2_exp=radix2_exp, window_type=WindowType.HANN,
                                     slide_length=(1 << radix2_exp) // 4,
                                     data_type=SpectralDataType.POWER,
                                     filter_bank_type=SpectralFilterBankType.LINEAR,
                                     filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                     filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newLinear']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self.num = self.get_band_num()
        self._is_created = True


class Mel(SpectrogramBase):
    """
    Simple mel spectrogram class.

    It can create simple mel spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Mel class fixed parameter:
        * fre range: 0 - 16000(samplate/2)
        * slide_length: ``fft_length // 4``
        * window_type: HANN
        * data_type: POWER
        * style_type: SLANEY
        * normal_type: NONE

    Parameters
    ----------
    num: int
        Number of mel bins. Usuall use 128/256/...
        
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::


        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get mel spectrogram data

        >>> from audioflux.spectrogram import Mel
        >>> spec_obj = Mel(num=128, samplate=32000, radix2_exp=12)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  y_coords=spec_obj.y_coords(),
        ...                  scale='log', title='Mel')
    """

    def __init__(self, num=128, samplate=32000, radix2_exp=12):
        super(Mel, self).__init__(num=num, samplate=samplate, low_fre=0.0, high_fre=16000.0,
                                  bin_per_octave=12,
                                  radix2_exp=radix2_exp, window_type=WindowType.HANN,
                                  slide_length=(1 << radix2_exp) // 4,
                                  data_type=SpectralDataType.POWER,
                                  filter_bank_type=SpectralFilterBankType.MEL,
                                  filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                  filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newMel']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class Bark(SpectrogramBase):
    """
    Simple bark spectrogram class.

    It can create simple bark spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Bark class fixed parameter:
        * fre range: 0 - 16000(samplate/2)
        * slide_length: ``fft_length // 4``
        * window_type: HANN
        * data_type: POWER
        * style_type: SLANEY
        * normal_type: NONE

    Parameters
    ----------
    num: int
        Number of bark bins. Usuall use 128/256/...
        
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::


        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get bark spectrogram data

        >>> from audioflux.spectrogram import Bark
        >>> spec_obj = Bark(num=128, samplate=32000, radix2_exp=12)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  y_coords=spec_obj.y_coords(),
        ...                  scale='log', title='Bark')
    """

    def __init__(self, num=128, samplate=32000, radix2_exp=12):
        super(Bark, self).__init__(num=num,
                                   samplate=samplate,
                                   low_fre=0.0,
                                   high_fre=16000.0,
                                   bin_per_octave=12,
                                   radix2_exp=radix2_exp,
                                   window_type=WindowType.HANN,
                                   slide_length=(1 << radix2_exp) // 4,
                                   data_type=SpectralDataType.POWER,
                                   filter_bank_type=SpectralFilterBankType.BARK,
                                   filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                   filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newBark']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class Erb(SpectrogramBase):
    """
    Simple erb spectrogram class.

    It can create simple erb spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Erb class fixed parameter:
        * fre range: 0 - 16000(samplate/2)
        * slide_length: ``fft_length // 4``
        * window_type: HANN
        * data_type: POWER
        * style_type: SLANEY
        * normal_type: NONE

    Parameters
    ----------
    num: int
        Number of erb bins. Usuall use 128/256/...
        
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::


        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get erb spectrogram data

        >>> from audioflux.spectrogram import Erb
        >>> spec_obj = Erb(num=128, samplate=32000, radix2_exp=12)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  y_coords=spec_obj.y_coords(),
        ...                  scale='log', title='Erb')
    """

    def __init__(self, num=128, samplate=32000, radix2_exp=12):
        super(Erb, self).__init__(num=num, samplate=samplate, low_fre=0.0, high_fre=16000.0,
                                  bin_per_octave=12,
                                  radix2_exp=radix2_exp, window_type=WindowType.HANN,
                                  slide_length=(1 << radix2_exp) // 4,
                                  data_type=SpectralDataType.POWER,
                                  filter_bank_type=SpectralFilterBankType.ERB,
                                  filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                  filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newErb']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class Chroma(SpectrogramBase):
    """
    Simple chroma spectrogram class.

    It can create simple chroma spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Chroma class fixed parameter:
        * chroma num: 12
        * fre range: 0 - 16000(samplate/2)
        * bin_per_octave: 12
        * slide_length: ``fft_length // 4``
        * window_type: HANN
        * data_type: POWER
        * style_type: SLANEY

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::

        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get erb spectrogram data

        >>> from audioflux.spectrogram import Chroma
        >>> spec_obj = Chroma()
        >>> spec_arr = spec_obj.spectrogram(audio_arr)

        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  title='Chroma')
    """

    def __init__(self, samplate=32000, radix2_exp=12):
        super(Chroma, self).__init__(num=12, samplate=samplate, low_fre=note_to_hz('C1'), high_fre=None,
                                     bin_per_octave=12,
                                     radix2_exp=radix2_exp, window_type=WindowType.HANN,
                                     slide_length=(1 << radix2_exp) // 4,
                                     data_type=SpectralDataType.POWER,
                                     filter_bank_type=SpectralFilterBankType.CHROMA,
                                     filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                     filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newChroma']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class Deep(SpectrogramBase):
    """
    Simple deep spectrogram class.

    It can create simple deep spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    Deep class fixed parameter:
        * fre range: 32.703(C1) - 16000(samplate/2)
        * slide_length: ``fft_length // 4``
        * window_type: HAMM
        * data_type: POWER
        * style_type: SLANEY

    Parameters
    ----------
    num: int:
        Number of deep bins.
        
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::

        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get deep spectrogram data

        >>> from audioflux.spectrogram import Deep
        >>> spec_obj = Deep(num=84)
        >>> spec_arr = spec_obj.spectrogram(audio_arr)

        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot(ncols=3, fig_height=10)
        >>> pt.add_spec_data(spec_arr[0],
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  title='Deep 0', col_idx=0)
        >>> pt.add_spec_data(spec_arr[1],
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  title='Deep 1', col_idx=1)
        >>> pt.add_spec_data(spec_arr[2],
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  title='Deep 2', col_idx=2)
    """

    def __init__(self, num, samplate=32000, radix2_exp=12):
        # WindowType 'hamm' ['hamm','hann','rect']
        # dataType 'power' mag||power
        # deepOrder 1 [1,2,3,4]
        # lowFre 32.703 highFre不起作用
        super(Deep, self).__init__(num=num, samplate=samplate, low_fre=note_to_hz('C1'), high_fre=16000.0,
                                   bin_per_octave=12,
                                   radix2_exp=radix2_exp, window_type=WindowType.HAMM,
                                   slide_length=(1 << radix2_exp) // 4,
                                   data_type=SpectralDataType.POWER,
                                   filter_bank_type=SpectralFilterBankType.DEEP,
                                   filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                   filter_normal_type=SpectralFilterBankNormalType.NONE)
        fn = self._lib['spectrogramObj_newDeep']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class DeepChroma(SpectrogramBase):
    """
    Simple DeepChroma spectrogram class.

    It can create simple DeepChroma spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``Spectrogram`` class to create.

    DeepChroma class fixed parameter:
        * chroma num: 12
        * low_fre: 32.703(C1)
        * slide_length: ``fft_length // 4``
        * window_type: HAMM
        * data_type: POWER
        * style_type: SLANEY

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio
        
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    See Also
    --------
    Spectrogram

    Examples
    --------
    .. plot::

        Get 220hz audio

        >>> import audioflux
        >>> path = audioflux.utils.sample_path('220')
        >>> audio_arr, sr = audioflux.read(path)

        Get DeepChroma spectrogram data

        >>> from audioflux.spectrogram import DeepChroma
        >>> spec_obj = DeepChroma()
        >>> spec_arr = spec_obj.spectrogram(audio_arr)


        Show plot

        >>> from audioflux.display import Plot
        >>> audio_len = audio_arr.shape[-1]
        >>> pt = Plot()
        >>> pt.add_spec_data(spec_arr,
        ...                  x_coords=spec_obj.x_coords(audio_len),
        ...                  title='DeepChroma')
    """

    def __init__(self, samplate=32000, radix2_exp=12):
        super(DeepChroma, self).__init__(num=12, samplate=samplate, low_fre=note_to_hz('C1'), high_fre=16000.0,
                                         bin_per_octave=12,
                                         radix2_exp=radix2_exp, window_type=WindowType.HAMM,
                                         slide_length=(1 << radix2_exp) // 4,
                                         data_type=SpectralDataType.POWER,
                                         filter_bank_type=SpectralFilterBankType.DEEP_CHROMA,
                                         filter_style_type=SpectralFilterBankStyleType.SLANEY,
                                         filter_normal_type=SpectralFilterBankNormalType.NONE)

        fn = self._lib['spectrogramObj_newDeepChroma']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectrogram)), c_int, c_int, POINTER(c_int)]
        fn(self._obj,
           c_int(self.samplate),
           c_int(self.radix2_exp),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True
