import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio

__all__ = ["HPSS"]


class OpaqueHPSS(Structure):
    _fields_ = []


class HPSS(Base):
    """
    HPSS - Median filtering, NMF algorithm.

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    h_order: int
        h order

    p_order: int
        p order

    Examples
    --------

    Get a chord audio file with a metronome

    >>> import audioflux as af
    >>> audio_arr, sr = af.read(af.utils.sample_path('chord_metronome1'))

    Create HPSS object and extrct h/p data

    >>> from audioflux.type import WindowType
    >>> radix2_exp = 11
    >>> slide_length = (1 << radix2_exp) // 4
    >>> hpss_obj = af.HPSS(radix2_exp=radix2_exp, window_type=WindowType.HAMM,
    >>>                    slide_length=slide_length, h_order=21, p_order=31)
    >>> h_arr, p_arr = hpss_obj.hpss(audio_arr)

    Disable Plot of Linear spectrogram

    >>> import numpy as np
    >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
    >>> bft_obj = af.BFT(num=2049, radix2_exp=12, samplate=sr,
    >>>                  scale_type=SpectralFilterBankScaleType.LINEAR,
    >>>                  data_type=SpectralDataType.POWER)
    >>> audio_arr = audio_arr[:len(h_arr)]
    >>> origin_spec_arr = bft_obj.bft(audio_arr, result_type=1)
    >>> h_spec_arr = bft_obj.bft(h_arr, result_type=1)
    >>> p_spec_arr = bft_obj.bft(p_arr, result_type=1)
    >>> origin_spec_arr = af.utils.power_to_abs_db(origin_spec_arr)
    >>> h_spec_arr = af.utils.power_to_abs_db(h_spec_arr)
    >>> p_spec_arr = af.utils.power_to_abs_db(p_spec_arr)

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[0]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(origin_spec_arr, axes=ax,
    >>>           x_coords=bft_obj.x_coords(audio_len),
    >>>           y_coords=bft_obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='Origin Linear Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(h_spec_arr, axes=ax,
    >>>           x_coords=bft_obj.x_coords(audio_len),
    >>>           y_coords=bft_obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='h_order Linear Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(p_spec_arr, axes=ax,
    >>>           x_coords=bft_obj.x_coords(audio_len),
    >>>           y_coords=bft_obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='p_order Linear Spectrogram')
    >>> fig.colorbar(img, ax=ax)

    """

    def __init__(self, radix2_exp=12, window_type=WindowType.HAMM,
                 slide_length=1024, h_order=21, p_order=31):
        super(HPSS, self).__init__(pointer(OpaqueHPSS()))

        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length
        self.h_order = h_order
        self.p_order = p_order

        fn = self._lib['hpssObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueHPSS)),
                       c_int,
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]

        fn(self._obj,
               c_int(self.radix2_exp),
               pointer(c_int(self.window_type.value)),
               pointer(c_int(self.slide_length)),
               pointer(c_int(self.h_order)),
               pointer(c_int(self.p_order)))
        self._is_created = True

    def cal_data_length(self, data_length):
        """
        Calculate the data length.

        Parameters
        ----------
        data_length: int
            Input array length

        Returns
        -------
        out: int
        """
        fn = self._lib['hpssObj_calDataLength']
        fn.argtypes = [POINTER(OpaqueHPSS), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def hpss(self, data_arr) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Compute the hpss

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]

        Returns
        -------
        h_arr: np.ndarray [shape=(n,)]
        p_arr: np.ndarray [shape=(n,)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)

        fn = self._lib['hpssObj_hpss']
        fn.argtypes = [POINTER(OpaqueHPSS),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        data_len = data_arr.size

        new_data_len = self.cal_data_length(data_len)
        h_arr = np.zeros(new_data_len, dtype=np.float32)
        p_arr = np.zeros(new_data_len, dtype=np.float32)

        fn(self._obj, data_arr, c_int(data_len), h_arr, p_arr)
        return h_arr, p_arr

    def debug(self):
        fn = self._lib['hpssObj_debug']
        fn.argtypes = [POINTER(OpaqueHPSS)]
        fn(self._obj)

    def __del__(self):
        if self._is_created:
            fn = self._lib['hpssObj_free']
            fn.argtypes = [POINTER(OpaqueHPSS)]
            fn.restype = c_void_p
            fn(self._obj)
