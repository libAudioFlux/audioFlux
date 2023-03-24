import numpy as np
import ctypes
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p, c_size_t
from audioflux.type import SpectralNoveltyMethodType, SpectralNoveltyDataType
from audioflux.base import Base
from audioflux.utils import ascontiguous_swapaxex, format_channel, revoke_channel

__all__ = ["Spectral"]


class OpaqueSpectral(Structure):
    _fields_ = []


class Spectral(Base):
    """
    Spectrum feature, supports all spectrum types.

    Parameters
    ----------
    num: int
        Number of frequency bins to generate. It must be the same as the
        `num` parameter of the transformation (same as the spectrogram matrix).

    fre_band_arr: np.ndarray [shape=(n_fre,)]
        The array of frequency bands. Obtained by calling the `get_fre_band_arr()` method of the transformation.
    """

    def __init__(self, num, fre_band_arr):
        super(Spectral, self).__init__(pointer(OpaqueSpectral()))

        self.num = num
        self.fre_band_arr = np.asarray(fre_band_arr, dtype=np.float32, order='C')

        self.time_length = 0
        fn = self._lib['spectralObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueSpectral)),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]
        fn(self._obj, c_int(self.num), self.fre_band_arr)
        self._is_created = True

    def set_time_length(self, time_length):
        """
        Set time length

        Parameters
        ----------
        time_length: int
        """

        fn = self._lib['spectralObj_setTimeLength']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            c_int
        ]
        fn(self._obj, c_int(time_length))
        self.time_length = time_length

    def set_edge(self, start, end):
        """
        Set edge

        Parameters
        ----------
        start: int
            0 ~ end
        end: int
            start ~ num-1
        """
        if not 0 <= start < end:
            raise ValueError(f'start={start} must be in range [0, {end})')
        if not start < end <= self.num - 1:
            raise ValueError(f'start={end} must be in range ({start}, {self.num - 1}]')

        fn = self._lib['spectralObj_setEdge']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            c_int,
            c_int,
        ]
        fn(self._obj, c_int(start), c_int(end))

    def set_edge_arr(self, index_arr):
        """
        Set edge array

        Parameters
        ----------
        index_arr: np.ndarray [shape=(n,), dtype=np.int32]
            fre index array
        """

        index_arr = np.asarray(index_arr, dtype=np.int32, order='C')
        if index_arr.ndim != 1:
            raise ValueError(f'index_arr must be a 1D array.')

        fn = self._lib['spectralObj_setEdgeArr']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
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

    def flatness(self, m_data_arr):
        """
        Compute the spectral flatness feature.

        :math:`\qquad flatness=\\frac{\\left ( \prod_{k=b_1}^{b_2} s_k  \\right)^{ \\frac{1}{b_2-b_1} } } {\\frac{1}{b_2-b_1} \sum_{ k=b_1 }^{b_2} s_k}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        flatness: np.ndarray [shape=(..., time)]
            flatness frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract flatness feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> flatness_arr = spectral_obj.flatness(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(flatness_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, flatness_arr, axes=ax[1], label='flatness')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_flatness']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * f
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def flux(self, m_data_arr, step=1, p=2, is_positive=False, is_no_exp=True, tp=0):
        """
        Compute the spectral flux feature.

        :math:`\qquad flux(t)=\\left( \sum_{k=b_1}^{b_2} |s_k(t)-s_k(t-1) |^{p}  \\right)^{\\frac{1}{p}}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        * In general :math:`s_k(t) \geq s_k(t-1)` participate in the calculation

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

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
        flux: np.ndarray [shape=(..., time)]
            flux frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract flux feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> flux_arr = spectral_obj.flux(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(flux_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, flux_arr, axes=ax[1], label='flux')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_flux']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int, c_float, c_int,
            POINTER(c_int), POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_time, n_fre = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_time, dtype=np.float32)
            fn(self._obj,
               m_data_arr,
               c_int(step),
               c_float(p),
               c_int(int(is_positive)),
               pointer(c_int(int(is_no_exp))),
               pointer(c_int(tp)),
               ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_time), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj,
                   m_data_arr[i],
                   c_int(step),
                   c_float(p),
                   c_int(int(is_positive)),
                   pointer(c_int(int(is_no_exp))),
                   pointer(c_int(tp)),
                   ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def rolloff(self, m_data_arr, threshold=0.95):
        """
        Compute the spectral rolloff feature.

        :math:`\qquad \sum_{k=b_1}^{i}|s_k| \geq \eta \sum_{k=b_1}^{b_2}s_k`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * :math:`\eta \in (0,1)`, generally take 0.95 or 0.85, satisfy the condition :math:`i` get :math:`f_i` rolloff frequency

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        threshold: float, [0,1]
            rolloff threshold. Generally take 0.95 or 0.85.

        Returns
        -------
        rolloff: np.ndarray [shape=(..., time)]
            rolloff frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract rolloff feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> rolloff_arr = spectral_obj.rolloff(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(rolloff_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, rolloff_arr, axes=ax[1], label='rolloff')
        """
        if not 0 <= threshold <= 1:
            raise ValueError(f'threshold must be 0 or 1')

        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_rolloff']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_float(threshold), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_float(threshold), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def centroid(self, m_data_arr):
        """
        Compute the spectral centroid feature.

        :math:`\qquad \mu_1=\\frac{\sum_{ k=b_1 }^{b_2} f_ks_k } {\sum_{k=b_1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`f_k` is in Hz
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        centroid: np.ndarray [shape=(..., time)]
            centroid frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract centroid feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> centroid_arr = spectral_obj.centroid(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(centroid_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, centroid_arr, axes=ax[1], label='centroid')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_centroid']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        spread: np.ndarray [shape=(..., time)]
            spread frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract spread feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> spread_arr = spectral_obj.spread(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(spread_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, spread_arr, axes=ax[1], label='spread')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_spread']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]

        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        skewness: np.ndarray [shape=(..., time)]
            skewness frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract skewness feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> skewness_arr = spectral_obj.skewness(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(skewness_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, skewness_arr, axes=ax[1], label='skewness')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_skewness']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]

        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        kurtosis: np.ndarray [shape=(..., time)]
            kurtosis frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract kurtosis feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> kurtosis_arr = spectral_obj.kurtosis(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(kurtosis_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, kurtosis_arr, axes=ax[1], label='kurtosis')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_kurtosis']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        is_norm: bool
            Whether to norm

        Returns
        -------
        entropy: np.ndarray [shape=(..., time)]
            entropy frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract entropy feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> entropy_arr = spectral_obj.entropy(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(entropy_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, entropy_arr, axes=ax[1], label='entropy')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_entropy']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(int(is_norm)), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(int(is_norm)), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def crest(self, m_data_arr):
        """
        Compute the spectral crest feature.

        :math:`\qquad crest =\\frac{max(s_{k\in_{[b_1,b_2]} }) } {\\frac{1}{b_2-b_1} \sum_{ k=b_1 }^{b_2} s_k}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        crest: np.ndarray [shape=(..., time)]
            crest frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract crest feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> crest_arr = spectral_obj.crest(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(crest_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, crest_arr, axes=ax[1], label='crest')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_crest']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        slope: np.ndarray [shape=(..., time)]
            slope frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract slope feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> slope_arr = spectral_obj.slope(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(slope_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, slope_arr, axes=ax[1], label='slope')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_slope']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def decrease(self, m_data_arr):
        """
        Compute the spectral decrease feature.

        :math:`\qquad decrease=\\frac { \sum_{k=b_1+1}^{b_2} \\frac {s_k-s_{b_1}}{k-1} } { \sum_{k=b_1+1}^{b_2} s_k }`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        decrease: np.ndarray [shape=(..., time)]
            decrease frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract decrease feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> decrease_arr = spectral_obj.decrease(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(decrease_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, decrease_arr, axes=ax[1], label='decrease')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_decrease']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        p: int, 1 or 2
            norm: 1 abs; 2 pow

        Returns
        -------
        band_width: np.ndarray [shape=(..., time)]
            band_width frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract band_width feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> band_width_arr = spectral_obj.band_width(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(band_width_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, band_width_arr, axes=ax[1], label='band_width')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_bandWidth']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_float(p), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_float(p), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def rms(self, m_data_arr):
        """
        Compute the spectral rms feature.

        :math:`\qquad rms=\sqrt{ \\frac{1}{N} \sum_{n=1}^N x^2[n] }=\sqrt {\\frac{1}{N^2}\sum_{m=1}^N |X[m]|^2 }`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        rms: np.ndarray [shape=(.... time)]
            rms frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract rms feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> rms_arr = spectral_obj.rms(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(rms_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, rms_arr, axes=ax[1], label='rms')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_rms']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def energy(self, m_data_arr, is_log=False, gamma=10.):
        """
        Compute the spectral energy feature.

        :math:`\qquad energy=\sum_{n=1}^N x^2[n] =\\frac{1}{N}\sum_{m=1}^N |X[m]|^2`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        is_log: bool
            Whether to log

        gamma: float
            energy gamma value.

        Returns
        -------
        energy: np.ndarray [shape=(..., time)]
            energy frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract energy feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> energy_arr = spectral_obj.energy(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(energy_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, energy_arr, axes=ax[1], label='energy')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_energy']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, int(is_log), gamma, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], int(is_log), gamma, ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def hfc(self, m_data_arr):
        """
        Compute the spectral hfc feature.

        :math:`\qquad hfc(t)=\\frac{\sum_{k=b_1}^{b_2} s_k(t)k }{b_2-b_1+1}`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        hfc: np.ndarray [shape=(..., time)]
            hfc frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract hfc feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> hfc_arr = spectral_obj.hfc(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(hfc_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, hfc_arr, axes=ax[1], label='hfc')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_hfc']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def sd(self, m_data_arr, step=1, is_positive=False):
        """
        Compute the spectral sd feature.

        :math:`\qquad sd(t)=flux(t)`

        satisfies the calculation of :math:`s_k(t) \ge s_k(t-1)`, :math:`p=2`，the result is not :math:`1/p`

        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * flux: `Spectral.flux`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        step: int
            Compute time axis steps, like 1/2/3/...

        is_positive: bool
            Whether to set negative numbers to 0

        Returns
        -------
        sd: np.ndarray [shape=(..., time)]
            sd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract sd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> sd_arr = spectral_obj.sd(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(sd_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, sd_arr, axes=ax[1], label='sd')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_sd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(step), c_int(int(is_positive)), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(step), c_int(int(is_positive)), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def sf(self, m_data_arr, step=1, is_positive=False):
        """
        Compute the spectral sf feature.

        :math:`\qquad sf(t)=flux(t)`

        satisfies the calculation of :math:`s_k(t) \ge s_k(t-1)`, :math:`p=1`

        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum
        * flux: `Spectral.flux`

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        step: int
            Compute time axis steps, like 1/2/3/...

        is_positive: bool
            Whether to set negative numbers to 0

        Returns
        -------
        sf: np.ndarray [shape=(..., time)]
            sf frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract sf feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> sf_arr = spectral_obj.sf(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(sf_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, sf_arr, axes=ax[1], label='sf')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_sf']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(step), c_int(int(is_positive)), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(step), c_int(int(is_positive)), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def mkl(self, m_data_arr, tp=0):
        """
        Compute the spectral mkl feature.

        :math:`\qquad mkl(t)=\sum_{k=b_1}^{b_2} \log\\left(1+ \cfrac {s_k(t)}{s_k(t-1)} \\right)`

        * :math:`b_1` and :math:`b_2`: the frequency band bin boundaries
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        tp: int, 0 or 1
            0 sum 1 mean

        Returns
        -------
        mkl: np.ndarray [shape=(..., time)]
            mkl frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract mkl feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> mkl_arr = spectral_obj.mkl(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(mkl_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, mkl_arr, axes=ax[1], label='mkl')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_mkl']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(tp), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(tp), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        m_phase_arr: np.ndarray [shape=(..., fre, time)]
            Phase data.

        Returns
        -------
        pd: np.ndarray [shape=(..., time)]
            pd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> phase_arr = af.utils.get_phase(spec_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract pd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> pd_arr = spectral_obj.pd(spec_arr, phase_arr)
        """

        if m_data_arr.shape != m_phase_arr.shape:
            raise ValueError(f'm_data_arr and m_phase_arr must be the same shape')

        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)

        fn = self._lib['spectralObj_pd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, m_phase_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            m_phase_arr, _ = format_channel(m_phase_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_phase_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        m_phase_arr: np.ndarray [shape=(..., fre, time)]
            Phase data.

        Returns
        -------
        wpd: np.ndarray [shape=(..., time)]
            wpd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> phase_arr = af.utils.get_phase(spec_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract wpd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> wpd_arr = spectral_obj.wpd(spec_arr, phase_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(wpd_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, wpd_arr, axes=ax[1], label='wpd')
        """

        if m_data_arr.shape != m_phase_arr.shape:
            raise ValueError(f'm_data_arr and m_phase_arr must be the same shape')
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)

        fn = self._lib['spectralObj_wpd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, m_phase_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            m_phase_arr, _ = format_channel(m_phase_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_phase_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def nwpd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral nwpd feature.

        :math:`\qquad nwpd(t)= \\frac {wpd} {\mu_s}`

        * wpd: `Spectral.wpd`
        * :math:`\mu_s`: the mean of :math:`s_k(t)`
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        m_phase_arr: np.ndarray [shape=(..., fre, time)]
            Phase data.

        Returns
        -------
        nwpd: np.ndarray [shape=(..., time)]
            nwpd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> phase_arr = af.utils.get_phase(spec_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract nwpd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> nwpd_arr = spectral_obj.nwpd(spec_arr, phase_arr)
        """

        if m_data_arr.shape != m_phase_arr.shape:
            raise ValueError(f'm_data_arr and m_phase_arr must be the same shape')
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)

        fn = self._lib['spectralObj_nwpd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, m_phase_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            m_phase_arr, _ = format_channel(m_phase_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_phase_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        m_phase_arr: np.ndarray [shape=(..., fre, time)]
            Phase data.

        Returns
        -------
        cd: np.ndarray [shape=(..., time)]
            cd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> phase_arr = af.utils.get_phase(spec_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract cd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> cd_arr = spectral_obj.cd(spec_arr, phase_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(cd_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, cd_arr, axes=ax[1], label='cd')
        """

        if m_data_arr.shape != m_phase_arr.shape:
            raise ValueError(f'm_data_arr and m_phase_arr must be the same shape')
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)

        fn = self._lib['spectralObj_cd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, m_phase_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            m_phase_arr, _ = format_channel(m_phase_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_phase_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def rcd(self, m_data_arr, m_phase_arr):
        """
        Compute the spectral rcd feature.

        :math:`\qquad rcd(t)=cd`

        participate in the sum calculation when :math:`s_k(t) \geq s_k(t-1)` is satisfied

        * cd: `Spectral.cd`
        * :math:`s_k`: the spectrum value, which can be magnitude spectrum or power spectrum

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        m_phase_arr: np.ndarray [shape=(..., fre, time)]
            Phase data.

        Returns
        -------
        rcd: np.ndarray [shape=(..., time)]
            rcd frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> phase_arr = af.utils.get_phase(spec_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract rcd feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> rcd_arr = spectral_obj.rcd(spec_arr, phase_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(rcd_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, rcd_arr, axes=ax[1], label='rcd')
        """

        if m_data_arr.shape != m_phase_arr.shape:
            raise ValueError(f'm_data_arr and m_phase_arr must be the same shape')
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_phase_arr = ascontiguous_swapaxex(m_phase_arr, -1, -2)

        fn = self._lib['spectralObj_rcd']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, m_phase_arr, ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            m_phase_arr, _ = format_channel(m_phase_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_phase_arr[i], ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def broadband(self, m_data_arr, threshold=0):
        """
        Compute the spectral broadband feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        threshold: float, [0,1]
            broadband threshold

        Returns
        -------
        broadband: np.ndarray [shape=(..., time)]
            broadband frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract broadband feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> broadband_arr = spectral_obj.broadband(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(broadband_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, broadband_arr, axes=ax[1], label='broadband')
        """
        if not 0 <= threshold <= 1:
            raise ValueError(f'threshold must be 0 or 1')

        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_broadband']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_float(threshold), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_float(threshold), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

        return ret_arr

    def novelty(self, m_data_arr, step=1, threshold=0.,
                method_type=SpectralNoveltyMethodType.SUB,
                data_type=SpectralNoveltyDataType.VALUE):
        """
        Compute the spectral novelty feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

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
        novelty: np.ndarray [shape=(..., time)]
            Novelty frequency per time step.

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract novelty feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> novelty_arr = spectral_obj.novelty(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(novelty_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, novelty_arr, axes=ax[1], label='novelty')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_novelty']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            POINTER(c_int),
            POINTER(c_int),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(step), c_float(threshold),
               pointer(c_int(method_type.value)),
               pointer(c_int(data_type.value)),
               ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(step), c_float(threshold),
                   pointer(c_int(method_type.value)),
                   pointer(c_int(data_type.value)),
                   ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        is_norm: bool
            Whether to norm

        Returns
        -------
        eef: np.ndarray [shape=(..., time)]
            eef frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract eef feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> eef_arr = spectral_obj.eef(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(eef_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, eef_arr, axes=ax[1], label='eef')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_eef']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(int(is_norm)), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(int(is_norm)), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)

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
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        is_norm: bool
            Whether to norm

        gamma: float
            Usually set is 1./10./20.etc, song is 0.5

        Returns
        -------
        eer: np.ndarray [shape=(..., time)]
            eer frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract eer feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> eer_arr = spectral_obj.eer(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(eer_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, eer_arr, axes=ax[1], label='eer')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_eer']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            c_int,
            c_float,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, c_int(int(is_norm)), c_float(gamma), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            ret_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], c_int(int(is_norm)), c_float(gamma), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def max(self, m_data_arr):
        """
        Compute the spectral max feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        val_arr: np.ndarray [shape=(..., time)]
            max value for each time

        fre_arr: np.ndarray [shape=(..., time)]
            max frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract max feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> max_val_arr, max_fre_arr = spectral_obj.max(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(max_val_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, max_val_arr, axes=ax[1], label='max_val')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_max']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            val_arr = np.zeros(n_len, dtype=np.float32)
            fre_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, val_arr, fre_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            val_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            fre_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], val_arr[i], fre_arr[i])
            val_arr = revoke_channel(val_arr, o_channel_shape, 1)
            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)
        return val_arr, fre_arr

    def mean(self, m_data_arr):
        """
        Compute the spectral mean feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        val_arr: np.ndarray [shape=(..., time)]
            mean value for each time

        fre_arr: np.ndarray [shape=(..., time)]
            mean frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract mean feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> mean_val_arr, mean_fre_arr = spectral_obj.mean(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(mean_val_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, mean_val_arr, axes=ax[1], label='mean_val')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_mean']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            val_arr = np.zeros(n_len, dtype=np.float32)
            fre_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, val_arr, fre_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            val_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            fre_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], val_arr[i], fre_arr[i])
            val_arr = revoke_channel(val_arr, o_channel_shape, 1)
            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)

        return val_arr, fre_arr

    def var(self, m_data_arr):
        """
        Compute the spectral var feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        val_arr: np.ndarray [shape=(..., time)]
            var value for each time

        fre_arr: np.ndarray [shape=(..., time)]
            var frequency for each time

        Examples
        --------

        Read chord audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('guitar_chord1')
        >>> audio_arr, sr = af.read(audio_path)

        Create BFT-Linear object and extract spectrogram

        >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
        >>> import numpy as np
        >>> bft_obj = af.BFT(num=2049, samplate=sr, radix2_exp=12, slide_length=1024,
        >>>                data_type=SpectralDataType.MAG,
        >>>                scale_type=SpectralFilterBankScaleType.LINEAR)
        >>> spec_arr = bft_obj.bft(audio_arr)
        >>> spec_arr = np.abs(spec_arr)

        Create Spectral object and extract var feature

        >>> spectral_obj = af.Spectral(num=bft_obj.num,
        >>>                            fre_band_arr=bft_obj.get_fre_band_arr())
        >>> n_time = spec_arr.shape[-1]  # Or use bft_obj.cal_time_length(audio_arr.shape[-1])
        >>> spectral_obj.set_time_length(n_time)
        >>> var_val_arr, var_fre_arr = spectral_obj.var(spec_arr)

        Display plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_plot, fill_wave
        >>> fig, ax = plt.subplots(nrows=2, sharex=True)
        >>> fill_wave(audio_arr, samplate=sr, axes=ax[0])
        >>> times = np.arange(0, len(var_val_arr)) * (bft_obj.slide_length / bft_obj.samplate)
        >>> fill_plot(times, var_val_arr, axes=ax[1], label='var_val')
        """
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['spectralObj_var']
        fn.argtypes = [
            POINTER(OpaqueSpectral),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        ]
        *_, n_len, m_len = m_data_arr.shape  # t * n
        if m_data_arr.ndim == 2:
            val_arr = np.zeros(n_len, dtype=np.float32)
            fre_arr = np.zeros(n_len, dtype=np.float32)
            fn(self._obj, m_data_arr, val_arr, fre_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            val_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            fre_arr = np.zeros((channel_num, n_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], val_arr[i], fre_arr[i])
            val_arr = revoke_channel(val_arr, o_channel_shape, 1)
            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)
        return val_arr, fre_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['spectralObj_free']
            fn.argtypes = [POINTER(OpaqueSpectral)]
            fn.restype = c_void_p
            fn(self._obj)
