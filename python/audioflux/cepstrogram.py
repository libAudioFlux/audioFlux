import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel, ascontiguous_swapaxex

__all__ = ['Cepstrogram']


class OpaqueCepstrogram(Structure):
    _fields_ = []


class Cepstrogram(Base):
    """
    Cepstrogram algorithm

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

     Examples
    --------

    Read guitar chord audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('guitar_chord2')
    >>> audio_arr, sr = af.read(audio_path)

    Extract Cepstrogram

    >>> from audioflux.type import ReassignType, WindowType
    >>> import numpy as np
    >>> obj = af.Cepstrogram(radix2_exp=12, samplate=sr)
    >>> cepstrums_arr, envelope_arr, details_arr = obj.cepstrogram(audio_arr)

    Show Cepstrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[-1]
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cepstrums_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Cepstrogram - Cepstrums')
    >>> fig.colorbar(img, ax=ax)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(envelope_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Cepstrogram - envelope')
    >>> fig.colorbar(img, ax=ax)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(details_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Cepstrogram - details')
    >>> fig.colorbar(img, ax=ax)
    >>>

    """

    def __init__(self, radix2_exp=12, samplate=32000, window_type=WindowType.RECT, slide_length=1024):
        super(Cepstrogram, self).__init__(pointer(OpaqueCepstrogram()))

        self.radix2_exp = radix2_exp
        self.slide_length = slide_length
        self.window_type = window_type
        self.samplate = samplate

        self.fft_length = 1 << radix2_exp

        fn = self._lib['cepstrogramObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueCepstrogram)),
                       c_int,
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           c_int(self.radix2_exp),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)))
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
        fn = self._lib['cepstrogramObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueCepstrogram), c_int]
        return fn(self._obj, c_int(data_length))

    def cepstrogram(self, data_arr, cep_num=4):
        """
        Get cepstrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data

        cep_num: int, 4~128
            formant estimate number

        Returns
        -------
        cepstrums: np.ndarray [shape=(..., fre, time), dtype=(np.float32)]
            The matrix of cepstrums

        envelope: np.ndarray [shape=(..., fre, time), dtype=(np.float32)]
            The matrix of envelope(formant)

        details: np.ndarray [shape=(..., fre, time), dtype=(np.float32)]
            The matrix of details(tone)
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        data_len = data_arr.shape[-1]

        fn = self._lib['cepstrogramObj_cepstrogram']
        fn.argtypes = [POINTER(OpaqueCepstrogram),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        time_len = self.cal_time_length(data_len)
        c_cep_num = c_int(cep_num)
        c_data_len = c_int(data_len)

        if data_arr.ndim == 1:
            size = (time_len, self.fft_length // 2 + 1)
            m_arr1 = np.zeros(size, dtype=np.float32)  # coef
            m_arr2 = np.zeros(size, dtype=np.float32)  # envelope
            m_arr3 = np.zeros(size, dtype=np.float32)  # tone
            fn(self._obj, c_cep_num, data_arr, c_data_len, m_arr1, m_arr2, m_arr3)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            size = (channel_num, time_len, self.fft_length // 2 + 1)
            m_arr1 = np.zeros(size, dtype=np.float32)  # envelope + tone
            m_arr2 = np.zeros(size, dtype=np.float32)  # envelope
            m_arr3 = np.zeros(size, dtype=np.float32)  # tone
            for i in range(channel_num):
                fn(self._obj, c_cep_num, data_arr[i], c_data_len,
                   m_arr1[i], m_arr2[i], m_arr3[i])
            m_arr1 = revoke_channel(m_arr1, o_channel_shape, 2)
            m_arr2 = revoke_channel(m_arr2, o_channel_shape, 2)
            m_arr3 = revoke_channel(m_arr3, o_channel_shape, 2)
        m_arr1 = ascontiguous_swapaxex(m_arr1, -2, -1)
        m_arr2 = ascontiguous_swapaxex(m_arr2, -2, -1)
        m_arr3 = ascontiguous_swapaxex(m_arr3, -2, -1)
        return m_arr1, m_arr2, m_arr3

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        y_coords = np.linspace(0, self.samplate / 2, int(self.fft_length / 2) + 2)
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
            free_fn = self._lib['cepstrogramObj_free']
            free_fn.argtypes = [POINTER(OpaqueCepstrogram)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
