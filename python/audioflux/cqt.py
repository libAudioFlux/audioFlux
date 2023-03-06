import warnings

import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.base import Base
from audioflux.type import (WindowType, SpectralFilterBankNormalType, SpectralDataType,
                            ChromaDataNormalType, CepstralRectifyType)
from audioflux.utils import check_audio, ascontiguous_T, note_to_hz

__all__ = [
    "CQT",
    "SimpleCQT"
]


class OpaqueCQT(Structure):
    _fields_ = []


class CQTBase(Base):
    def __init__(self, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
                 factor=1., beta=0., thresh=0.01,
                 window_type=WindowType.HANN, slide_length=None, is_continue=False,
                 normal_type=SpectralFilterBankNormalType.AREA,
                 is_scale=True):
        super(CQTBase, self).__init__(pointer(OpaqueCQT()))

        if bin_per_octave not in (12, 24, 36):
            raise ValueError(f'bin_per_octave={bin_per_octave} must be 12, 24 or 36')
        if num % bin_per_octave != 0:
            raise ValueError(f'num={num} must be an integer multiple of bin_per_octave={bin_per_octave}')

        self.num = num
        self.samplate = samplate
        self.low_fre = low_fre
        self.bin_per_octave = bin_per_octave
        self.factor = factor
        self.beta = beta
        self.thresh = thresh
        self.window_type = window_type
        self.slide_length = slide_length
        self.is_continue = is_continue
        self.normal_type = normal_type
        self.is_scale = is_scale

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        ``(data_length - fft_length) / slide_length + 1``

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['cqtObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueCQT), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def get_fft_length(self):
        """
        Get fft_length

        Returns
        -------
        out: int
        """
        fn = self._lib['cqtObj_getFFTLength']
        fn.argtypes = [POINTER(OpaqueCQT)]
        fn.restype = c_int
        return fn(self._obj)

    def get_fre_band_arr(self):
        """
        Get an array of frequency bands of CQT scales.

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        fn = self._lib['cqtObj_getFreBandArr']
        fn.argtypes = [POINTER(OpaqueCQT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_float * self.num).from_address(p), np.float32).copy()
        return ret

    def set_scale(self, flag):
        """
        Set scale

        Parameters
        ----------
        flag: int
        """
        fn = self._lib['cqtObj_setScale']
        fn.argtypes = [POINTER(OpaqueCQT), c_int]
        fn(self._obj, c_int(flag))

    def cqt(self, data_arr):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n)]
            Input audio data

        Returns
        -------
        out: np.ndarray [shape=(fre, time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)

        cqt_fn = self._lib['cqtObj_cqt']
        cqt_fn.argtypes = [POINTER(OpaqueCQT),
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                           c_int,
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                           np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                           ]

        data_len = data_arr.shape[0]
        time_len = self.cal_time_length(data_len)

        m_real_arr = np.zeros((time_len, self.num), dtype=np.float32)
        m_image_arr = np.zeros((time_len, self.num), dtype=np.float32)

        cqt_fn(self._obj, data_arr, c_int(data_len), m_real_arr, m_image_arr)

        ret_arr = m_real_arr + m_image_arr * 1j
        ret_arr = ascontiguous_T(ret_arr)
        return ret_arr

    def chroma(self, m_cqt_data, chroma_num=12,
               data_type=SpectralDataType.POWER,
               norm_type=ChromaDataNormalType.MAX):
        """
        Calculate the chroma matrix of CQT

        Parameters
        ----------
        m_cqt_data: np.ndarray [shape=(fre, time), dtype=np.complex]
            CQT spectrogram matrix, call the **cqt** method to get.

        chroma_num: int
            Number of chroma bins to produce

        data_type: SpectralDataType
            Data type of CQT spectrogram matrix

            See: `type.SpectralDataType`

        norm_type: ChromaDataNormalType
            Normalization type of chroma

        Returns
        -------
        out: np.ndarray [shape=(chroma_num, time)]
        """
        if not np.iscomplexobj(m_cqt_data):
            raise ValueError(f'm_cqt_data must be of type np.complex')
        if m_cqt_data.ndim != 2:
            raise ValueError(f'm_cqt_data is only defined for 2D arrays')

        m_cqt_data = ascontiguous_T(m_cqt_data)

        chroma_fn = self._lib['cqtObj_chroma']
        chroma_fn.argtypes = [POINTER(OpaqueCQT),
                              POINTER(c_int),
                              POINTER(c_int),
                              POINTER(c_int),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              ]

        t_len, _ = m_cqt_data.shape  # (time, fre)
        m_ret_arr = np.zeros((t_len, chroma_num), dtype=np.float32)

        m_real_arr = m_cqt_data.real.astype(dtype=np.float32)
        m_imag_arr = m_cqt_data.imag.astype(dtype=np.float32)

        chroma_fn(self._obj,
                  pointer(c_int(chroma_num)),
                  pointer(c_int(data_type.value)),
                  pointer(c_int(norm_type.value)),
                  m_real_arr, m_imag_arr, m_ret_arr)
        return ascontiguous_T(m_ret_arr)

    def cqcc(self, m_data_arr, cc_num, rectify_type) -> np.ndarray:
        """
        Compute the spectral cqcc feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            CQT spectrogram matrix, call the **cqt** method to get.

            - data_type:
                1. mag: ``np.abs(D)``
                2. power: ``np.abs(D) ** 2``
            - If data of `np.complex` type is passed in, the default is `power`

        cc_num: int
            Number of cc to produce

        rectify_type: CepstralRectifyType
            Rectify type

        Returns
        -------
        out: np.ndarray [shape=(cc_num, time)]
        """
        if m_data_arr.ndim != 2:
            raise ValueError(f'm_data_arr is only defined for 2D arrays')
        if np.iscomplexobj(m_data_arr):
            m_data_arr = np.abs(m_data_arr) ** 2
        m_data_arr = ascontiguous_T(m_data_arr)

        cqcc_fn = self._lib['cqtObj_cqcc']
        cqcc_fn.argtypes = [POINTER(OpaqueCQT),
                            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                            c_int,
                            POINTER(c_int),
                            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                            ]

        time_len, _ = m_data_arr.shape  # (time, fre)
        cqcc_arr = np.zeros((time_len, cc_num), dtype=np.float32)

        cqcc_fn(self._obj, m_data_arr, c_int(cc_num), pointer(c_int(rectify_type.value)), cqcc_arr)
        return ascontiguous_T(cqcc_arr)

    def cqhc(self, m_data_arr, hc_num=20) -> np.ndarray:
        """
        Compute the spectral cqhc feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            CQT spectrogram matrix, call the **cqt** method to get.

            - data_type:
                1. mag: ``np.abs(D)``
                2. power: ``np.abs(D) ** 2``
            - If data of `np.complex` type is passed in, the default is `power`

        hc_num: int
            Number of hc to produce

        Returns
        -------
        out: np.ndarray [shape=(hc_num, time)]
        """
        if m_data_arr.ndim != 2:
            raise ValueError(f'm_data_arr is only defined for 2D arrays')
        if np.iscomplexobj(m_data_arr):
            m_data_arr = np.abs(m_data_arr) ** 2
        m_data_arr = ascontiguous_T(m_data_arr)

        cqhc_fn = self._lib['cqtObj_cqhc']
        cqhc_fn.argtypes = [POINTER(OpaqueCQT),
                            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                            c_int,
                            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                            ]

        n_len, m_len = m_data_arr.shape
        m_form_arr = np.zeros((n_len, hc_num), dtype=np.float32)
        cqhc_fn(self._obj, m_data_arr, c_int(hc_num), m_form_arr)
        return ascontiguous_T(m_form_arr)

    def deconv(self, m_data_arr) -> (np.ndarray, np.ndarray):
        """
        Compute the spectral deconv feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(fre, time)]
            CQT spectrogram matrix, call the **cqt** method to get.

            - data_type:
                1. mag: ``np.abs(D)``
                2. power: ``np.abs(D) ** 2``
            - If data of `np.complex` type is passed in, the default is `mag`

        Returns
        -------
        m_tone_arr: np.ndarray [shape=(..., time)]
            The matrix of tone
        m_pitch_arr: np.ndarray [shape=(..., time)]
            The matrix of pitch
        """
        if m_data_arr.ndim != 2:
            raise ValueError(f'm_data_arr is only defined for 2D arrays')
        if np.iscomplexobj(m_data_arr):
            m_data_arr = np.abs(m_data_arr)
        m_data_arr = ascontiguous_T(m_data_arr)

        deconv_fn = self._lib['cqtObj_deconv']
        deconv_fn.argtypes = [POINTER(OpaqueCQT),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                              ]

        m_tone_arr = np.zeros_like(m_data_arr, dtype=np.float32)
        m_pitch_arr = np.zeros_like(m_data_arr, dtype=np.float32)
        deconv_fn(self._obj, m_data_arr, m_tone_arr, m_pitch_arr)

        m_tone_arr = ascontiguous_T(m_tone_arr)
        m_pitch_arr = ascontiguous_T(m_pitch_arr)
        return m_tone_arr, m_pitch_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate of CQT

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
        x_coords = np.linspace(0, data_length / self.samplate,
                               self.cal_time_length(data_length) + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['cqtObj_free']
            free_fn.argtypes = [POINTER(OpaqueCQT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)


class SimpleCQT(CQTBase):
    """
    Simple CQT spectrogram class.

    It can create simple CQT spectrogram, and only set a few basic parameters.
    If you want more parameter settings, use ``CQT`` class to create.

    SimpleCQT class fixed parameter:
        * bin_per_octave: 12
        * low_fre: 32.703(C1)
        * factor: 1
        * beta: 0
        * thresh: 0.01
        * window_type: HANN
        * slide_length: ``fft_length / 4``
        * normal_type: AREA
        * is_scale: 1

    Parameters
    ----------
    num: int
        Number of frequency bins to generate, starting at `low_fre`.

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)
    array([-5.5879354e-09, -9.3132257e-09,  0.0000000e+00, ...,
           3.2826858e-03,  3.2447521e-03,  3.0795704e-03], dtype=float32)

    Create SimpleCQT object

    >>> from audioflux.type import SpectralFilterBankNormalType
    >>> from audioflux.utils import note_to_hz
    >>> obj = af.SimpleCQT(num=84, samplate=sr, low_fre=note_to_hz('C1'))

    Extract CQT spectrogram

    >>> import numpy as np
    >>> spec_arr = obj.cqt(audio_arr)
    >>> spec_mag_arr = np.abs(spec_arr)
    array([[4.13306177e-01, 4.24238801e-01, 4.03985798e-01, ...,
            9.09681246e-03, 7.72366347e-03, 6.51519699e-03],
           [3.45392436e-01, 3.47385347e-01, 3.18555593e-01, ...,
            3.74186062e-03, 3.57810827e-03, 5.09629818e-03],
           [2.74609476e-01, 2.69073159e-01, 2.38167673e-01, ...,
            1.73897278e-02, 1.37754036e-02, 1.05445199e-02],
           ...,
           [6.24795386e-04, 4.99437414e-02, 1.60291921e-02, ...,
            1.42244404e-04, 1.27240157e-04, 7.22698495e-03],
           [8.64346686e-04, 3.18217799e-02, 2.55288873e-02, ...,
            1.96875044e-04, 1.22444399e-04, 7.08540343e-03],
           [5.32156206e-04, 3.65380235e-02, 1.94845423e-02, ...,
            1.01396305e-04, 4.90328348e-05, 6.83017401e-03]], dtype=float32)

    Show CQT spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[0]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_mag_arr, axes=ax,
    >>>           x_coords=obj.x_coords(audio_len),
    >>>           y_coords=obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='CQT Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    Extract Chroma-cqt data

    >>> chroma_arr = obj.chroma(spec_arr, chroma_num=12)
    array([[2.28795856e-01, 1.27659831e-02, 5.76458964e-03, ...,
            1.08233641e-03, 1.06033171e-03, 8.46851245e-02],
           ...,
           [3.88875484e-01, 4.20701923e-03, 1.00626342e-03, ...,
            7.09769898e-04, 3.13778268e-03, 2.85723597e-01]], dtype=float32)

    Show Chroma-CQT spectrogram plot

    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(chroma_arr, axes=ax,
    >>>           x_coords=obj.x_coords(audio_len),
    >>>           x_axis='time', y_axis='chroma',
    >>>           title='Chroma-CQT Spectrogram')
    >>> fig.colorbar(img, ax=ax)


    """

    def __init__(self, num=84, samplate=32000, low_fre=note_to_hz('C1')):
        super(SimpleCQT, self).__init__(num=num,
                                        samplate=samplate,
                                        low_fre=low_fre,
                                        bin_per_octave=12,
                                        factor=1., beta=0., thresh=0.01,
                                        window_type=WindowType.HANN, slide_length=None,
                                        normal_type=SpectralFilterBankNormalType.AREA,
                                        is_scale=True)

        fn = self._lib['cqtObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueCQT)),
                       c_int,
                       c_int,
                       c_float,
                       POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.samplate),
           c_float(self.low_fre),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True


class CQT(CQTBase):
    """
    Constant-Q transform (CQT)

    Parameters
    ----------
    num: int
        Number of frequency bins to generate, starting at `low_fre`.

        Usually: ``num = octave * bin_per_octave``, `default: 84 (7 * 12)`

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    bin_per_octave: int
        Number of bins per octave.

    factor: float
        Factor value

    beta: float
        Beta value

    thresh: float
        Thresh value

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    is_scale: bool
        Whether to use scale.

    See Also
    --------
    ST
    FST
    DWT
    WPT
    SWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Create CQT object

    >>> from audioflux.type import SpectralFilterBankNormalType
    >>> from audioflux.utils import note_to_hz
    >>> obj = af.CQT(num=84, samplate=sr, low_fre=note_to_hz('C1'), bin_per_octave=12,
    >>>              slide_length=1024, normal_type=SpectralFilterBankNormalType.AREA)

    Extract CQT spectrogram

    >>> import numpy as np
    >>> spec_arr = obj.cqt(audio_arr)
    >>> spec_mag_arr = np.abs(spec_arr)

    Show CQT spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[0]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_mag_arr, axes=ax,
    >>>           x_coords=obj.x_coords(audio_len),
    >>>           y_coords=obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='CQT Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    Extract Chroma-cqt data

    >>> chroma_arr = obj.chroma(spec_arr, chroma_num=12)
    array([[2.28795856e-01, 1.27659831e-02, 5.76458964e-03, ...,
            1.08233641e-03, 1.06033171e-03, 8.46851245e-02],
           ...,
           [3.88875484e-01, 4.20701923e-03, 1.00626342e-03, ...,
            7.09769898e-04, 3.13778268e-03, 2.85723597e-01]], dtype=float32)

    Show Chroma-CQT spectrogram plot

    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(chroma_arr, axes=ax,
    >>>           x_coords=obj.x_coords(audio_len),
    >>>           x_axis='time', y_axis='chroma',
    >>>           title='Chroma-CQT Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    """

    def __init__(self, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
                 factor=1., beta=0., thresh=0.01,
                 window_type=WindowType.HANN, slide_length=None,
                 normal_type=SpectralFilterBankNormalType.AREA,
                 is_scale=True):
        super(CQT, self).__init__(num=num, samplate=samplate, low_fre=low_fre,
                                  bin_per_octave=bin_per_octave,
                                  factor=factor, beta=beta, thresh=thresh,
                                  window_type=window_type, slide_length=slide_length,
                                  normal_type=normal_type,
                                  is_scale=is_scale)
        fn = self._lib['cqtObj_newWith']
        fn.argtypes = [POINTER(POINTER(OpaqueCQT)),
                       c_int,
                       POINTER(c_int),
                       POINTER(c_float),
                       POINTER(c_int),
                       POINTER(c_float),
                       POINTER(c_float),
                       POINTER(c_float),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           c_int(num),
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_int(self.bin_per_octave)),
           pointer(c_float(self.factor)),
           pointer(c_float(self.beta)),
           pointer(c_float(self.thresh)),
           pointer(c_int(self.window_type.value)),
           None if self.slide_length is None else pointer(c_int(self.slide_length)),
           pointer(c_int(int(self.is_continue))),
           pointer(c_int(self.normal_type.value)),
           pointer(c_int(int(self.is_scale))))
        self._is_created = True
