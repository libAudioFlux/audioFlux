import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.base import Base
from audioflux.type import NoveltyType
from audioflux.utils import ascontiguous_swapaxex, format_channel, revoke_channel

__all__ = ["Onset", "NoveltyParam"]


class OpaqueOnset(Structure):
    _fields_ = []


class NoveltyParam(Structure):
    """
    Onset's novelty parameters

    The parameters of the `Onset.novelty_type` corresponding method.

    Which parameter is used, please see: `Spectral`
    """
    _fields_ = [
        ("step", c_int),
        ("p", c_float),
        ("isPostive", c_int),
        ("isExp", c_int),
        ("type", c_int),
        ("threshold", c_float),
        ("isNorm", c_int),
        ("gamma", c_float),
    ]


class Onset(Base):
    """
    Onset - Spectrum flux, novelty, etc algorithm

    Parameters
    ----------
    time_length: int
        The length of the time axis of the Spectrogram matrix

    fre_length: int
        The length of the frequency axis of the Spectrogram matrix

    slide_length: int
        Sliding length, needs to be the same as Spectrogram

    samplate: int
        Sampling rate, needs to be the same as Spectrogram

    filter_order: int
        Filter order

    novelty_type: NoveltyType
        Novelty type

    Examples
    --------
    >>> import audioflux as af
    >>> audio_arr, sr = af.read(af.utils.sample_path('guitar_chord1'))

    >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
    >>> import numpy as np
    >>> bft_obj = af.BFT(num=128, samplate=sr, radix2_exp=12, slide_length=2048,
    >>>                  scale_type=SpectralFilterBankScaleType.MEL,
    >>>                  data_type=SpectralDataType.POWER)
    >>> spec_arr = bft_obj.bft(audio_arr)
    >>> spec_dB_arr = af.utils.power_to_db(np.abs(spec_arr))

    >>> from audioflux.type import NoveltyType
    >>> n_fre, n_time = spec_dB_arr.shape
    >>> onset_obj = af.Onset(time_length=n_time, fre_length=n_fre,
    >>>                      slide_length=bft_obj.slide_length, samplate=bft_obj.samplate,
    >>>                      novelty_type=NoveltyType.FLUX)
    >>> params = af.NoveltyParam(1, 2, 0, 1, 0, 0, 0, 1)
    >>> point_arr, evn_arr, time_arr, value_arr = onset_obj.onset(spec_dB_arr, novelty_param=params)

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec, fill_wave, fill_plot
    >>> audio_len = audio_arr.shape[-1]
    >>> fig, axes = plt.subplots(nrows=3, sharex=True)
    >>> img = fill_spec(spec_dB_arr, axes=axes[0],
    >>>                 x_coords=bft_obj.x_coords(audio_len),
    >>>                 y_coords=bft_obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Onset')
    >>>
    >>> ax = fill_wave(audio_arr, samplate=sr, axes=axes[1])
    >>>
    >>> times = np.arange(0, len(evn_arr)) * (bft_obj.slide_length / sr)
    >>> ax = fill_plot(times, evn_arr, axes=axes[2], label='Onset strength')
    >>> ax.vlines(time_arr, evn_arr.min(), evn_arr.max(), color='r', alpha=0.9,
    >>>           linestyle='--', label='Onsets')
    """

    def __init__(self, time_length, fre_length, slide_length, samplate=32000,
                 filter_order=1, novelty_type=NoveltyType.FLUX):
        super(Onset, self).__init__(pointer(OpaqueOnset()))

        self.time_length = time_length
        self.fre_length = fre_length
        self.samplate = samplate
        self.slide_length = slide_length
        self.filter_order = filter_order
        self.novelty_type = novelty_type

        fn = self._lib['onsetObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueOnset)), c_int, c_int, c_int,
                       POINTER(c_int), POINTER(c_int), POINTER(c_int)]

        fn(self._obj,
           c_int(self.time_length),
           c_int(self.fre_length),
           c_int(self.slide_length),
           pointer(c_int(self.samplate)),
           pointer(c_int(self.filter_order)),
           pointer(c_int(self.novelty_type.value)))
        self._is_created = True

    def onset(self, m_data_arr1, m_data_arr2=None, novelty_param=None, index_arr=None):
        """
        Compute onset

        Parameters
        ----------
        m_data_arr1: np.ndarray [shape=(..., fre, time)]
            Input spec data.

        m_data_arr2: np.ndarray [shape=(..., fre, time)] or None
            Input phase data. Provided when `novelty_type` is PD/WPD/NWPD/CD/RCD

        novelty_param: NoveltyParam or None
            The parameters of the `novelty_type` corresponding method.

            See: `NoveltyParam`

        index_arr: np.ndarray [shape=()] or None
            The index of frequency array

        Returns
        -------
        point_arr: np.ndarray [shape=(..., time)]
        evn_arr: np.ndarray [shape=(..., time)]
        time_arr: np.ndarray [shape=(..., time)]
        value_arr: np.ndarray [shape=(..., time)]
        """
        m_data_arr1 = np.asarray(m_data_arr1, dtype=np.float32, order='C')
        m_data_arr1 = ascontiguous_swapaxex(m_data_arr1, -1, -2)
        if m_data_arr2 is not None:
            m_data_arr2 = np.asarray(m_data_arr2, dtype=np.float32, order='C')
            m_data_arr2 = ascontiguous_swapaxex(m_data_arr2, -1, -2)
            if m_data_arr1.shape != m_data_arr2.shape:
                raise ValueError(f'm_data_arr1 and m_data_arr2 must be the same shape')

        if novelty_param is None:
            novelty_param = NoveltyParam(c_int(1), c_float(1), c_int(1), c_int(0),
                                         c_int(1), c_float(0), c_int(1), c_float(1))
        elif not isinstance(novelty_param, NoveltyParam):
            raise ValueError(f'novelty_param must be type of NoveltyParam')

        fn = self._lib['onsetObj_onset']
        fn.argtypes = [POINTER(OpaqueOnset),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       POINTER(c_void_p) if m_data_arr2 is None else np.ctypeslib.ndpointer(dtype=np.float32,
                                                                                            ndim=2,
                                                                                            flags='C_CONTIGUOUS'),
                       POINTER(NoveltyParam),
                       POINTER(c_void_p) if index_arr is None else np.ctypeslib.ndpointer(dtype=np.int32,
                                                                                          ndim=1,
                                                                                          flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        index_arr = None if index_arr is None else index_arr.astype(dtype=np.int32)
        index_length = 0 if index_arr is None else len(index_arr)

        if m_data_arr1.ndim == 2:
            n_time, n_fre = m_data_arr1.shape
            evn_arr = np.zeros(n_time, dtype=np.float32)
            point_arr = np.zeros(n_time, dtype=np.int32)
            point_len = fn(self._obj, m_data_arr1, m_data_arr2, novelty_param,
                           index_arr, c_int(index_length),
                           evn_arr, point_arr)
            point_arr = point_arr[:point_len]
            value_arr = evn_arr[point_arr]
        else:
            m_data_arr1, o_channel_shape = format_channel(m_data_arr1, 2)
            channel_num, n_time, n_fre = m_data_arr1.shape
            if m_data_arr2 is not None:
                m_data_arr2, _ = format_channel(m_data_arr2, 2)

            evn_arr = np.zeros((channel_num, n_time), dtype=np.float32)
            point_arr = np.zeros((channel_num, n_time), dtype=np.int32)
            value_arr = np.zeros((channel_num, n_time), dtype=np.float32)
            point_len = 0
            for i in range(channel_num):
                _m_data_arr2 = None if m_data_arr2 is None else m_data_arr2[i]
                _point_len = fn(self._obj, m_data_arr1[i], _m_data_arr2, novelty_param,
                                index_arr, c_int(index_length),
                                evn_arr[i], point_arr[i])
                point_len = max(point_len, _point_len)
                value_arr[i] = evn_arr[i, point_arr[i, :_point_len]]
            point_arr = point_arr[..., :point_len]
            value_arr = value_arr[..., :point_len]
            evn_arr = revoke_channel(evn_arr, o_channel_shape, 1)
            point_arr = revoke_channel(point_arr, o_channel_shape, 1)
            value_arr = revoke_channel(value_arr, o_channel_shape, 1)

        time_arr = 1.0 * point_arr * self.slide_length / self.samplate

        return point_arr, evn_arr, time_arr, value_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['onsetObj_free']
            fn.argtypes = [POINTER(OpaqueOnset)]
            fn.restype = c_void_p
            fn(self._obj)
