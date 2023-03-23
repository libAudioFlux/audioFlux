import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.base import Base
from audioflux.utils import ascontiguous_swapaxex, format_channel, revoke_channel

__all__ = ["Deconv"]


class OpaqueDeconv(Structure):
    _fields_ = []


class Deconv(Base):
    """
    Deconvolution for spectrum, supports all spectrum types.

    Parameters
    ----------
    num: int
        Number of frequency bins to generate. It must be the same as the
        `num` parameter of the transformation (same as the spectrogram matrix).

    Examples
    --------

    Get a 880Hz's audio file

    >>> import audioflux as af
    >>> sample_path = af.utils.sample_path('880')
    >>> audio_arr, sr = af.read(sample_path)

    Create BFT object and extract mel spectrogram

    >>> import numpy as np
    >>> from audioflux.type import SpectralFilterBankScaleType, SpectralDataType
    >>> bft_obj = af.BFT(num=128, radix2_exp=12, samplate=sr,
    >>>                  scale_type=SpectralFilterBankScaleType.MEL,
    >>>                  data_type=SpectralDataType.POWER)
    >>> spec_arr = bft_obj.bft(audio_arr)
    >>> spec_arr = np.abs(spec_arr)

    Create Deconv object and extract deconv

    >>> deconv_obj = af.Deconv(bft_obj.num)
    >>> deconv_obj.set_time_length(time_length=spec_arr.shape[1])
    >>> tone_arr, pitch_arr = deconv_obj.deconv(spec_arr)

    Display Deconv

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[-1]
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(tone_arr, axes=ax,
    >>>           x_coords=bft_obj.x_coords(audio_len), x_axis='time',
    >>>           title='Deconv Tone')
    >>> fig.colorbar(img, ax=ax)
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(pitch_arr, axes=ax,
    >>>           x_coords=bft_obj.x_coords(audio_len), x_axis='time',
    >>>           title='Deconv Pitch')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, num):
        super(Deconv, self).__init__(pointer(OpaqueDeconv()))

        self.num = num

        self.time_length = 0

        fn = self._lib['deconvObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueDeconv)), c_int]
        fn(self._obj, c_int(self.num))
        self._is_created = True

    def set_time_length(self, time_length):
        """
        Set time length

        Parameters
        ----------
        time_length: int
        """

        fn = self._lib['deconvObj_setTimeLength']
        fn.argtypes = [POINTER(OpaqueDeconv), c_int]
        fn(self._obj, c_int(time_length))
        self.time_length = time_length

    def deconv(self, m_data_arr):
        """
        Compute the spectral deconv feature.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time)]
            Spectrogram data.

        Returns
        -------
        m_tone_arr: np.ndarray [shape=(..., time)]
            The matrix of tone

        m_pitch_arr: np.ndarray [shape=(..., time)]
            The matrix of pitch
        """
        if m_data_arr.ndim < 2:
            raise ValueError(f"m_data_arr's dimensions must be greater than 1")
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)

        fn = self._lib['deconvObj_deconv']
        fn.argtypes = [POINTER(OpaqueDeconv),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        if m_data_arr.ndim == 2:
            m_tone_arr = np.zeros_like(m_data_arr, dtype=np.float32)
            m_pitch_arr = np.zeros_like(m_data_arr, dtype=np.float32)
            fn(self._obj, m_data_arr, m_tone_arr, m_pitch_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num, *_ = m_data_arr.shape

            m_tone_arr = np.zeros_like(m_data_arr, dtype=np.float32)
            m_pitch_arr = np.zeros_like(m_data_arr, dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_data_arr[i], m_tone_arr[i], m_pitch_arr[i])
            m_tone_arr = revoke_channel(m_tone_arr, o_channel_shape, 2)
            m_pitch_arr = revoke_channel(m_pitch_arr, o_channel_shape, 2)

        m_tone_arr = ascontiguous_swapaxex(m_tone_arr, -1, -2)
        m_pitch_arr = ascontiguous_swapaxex(m_pitch_arr, -1, -2)
        return m_tone_arr, m_pitch_arr

    def __del__(self):
        if self._is_created:
            fn = self._lib['deconvObj_free']
            fn.argtypes = [POINTER(OpaqueDeconv)]
            fn.restype = c_void_p
            fn(self._obj)
