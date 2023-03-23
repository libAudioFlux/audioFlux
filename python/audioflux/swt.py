import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.type import WaveletDiscreteType
from audioflux.base import Base
from audioflux.utils import check_audio, format_channel, revoke_channel

__all__ = ["SWT"]


class OpaqueSWT(Structure):
    _fields_ = []


class SWT(Base):
    """
    Stationary Wavelet Transform (SWT)

    Parameters
    ----------
    num: int
        Number of frequency bins to generate.

    fft_length: int
        fft length

    wavelet_type: WaveletDiscreteType
        Wavelet discrete type

        See: `type.WaveletDiscreteType`

        .. note::
            t1/t2 settings for wavelet_type:

            - DB: t1
                - 2~10/20/30/40
            - SYM: t1
                - 2~10/20/30
            - COIF: t1
                - 1/2/3/4/5
            - FK: t1
                - 4/6/8/14/18/22
            - BIOR/DMEY: t1.t2
                - 1.1/1.3/1.5
                - 2.2/2.4/2.6/2.8
                - 3.1/3.3/3.5/3.7/3.9
                - 4.4/5.5/6.8

    t1: int
        t1 value

    t2: int
        t2 value

    See Also
    --------
    CQT
    ST
    FST
    DWT
    WPT

    Examples
    --------
    .. plot::

        Read 220Hz audio data

        >>> import audioflux as af
        >>> audio_path = af.utils.sample_path('220')
        >>> audio_arr, sr = af.read(audio_path)
        >>> audio_len = 4096 * 5
        >>> audio_arr = audio_arr[:audio_len]

        Create SWT object

        >>> from audioflux.type import WaveletDiscreteType
        >>> obj = af.SWT(num=5, fft_length=audio_len,
        >>>              wavelet_type=WaveletDiscreteType.DB, t1=4, t2=0)

        Get SWT data

        >>> data_arr1, data_arr2 = obj.swt(audio_arr)

        Save SWT data

        >>> save_data_arr1 = data_arr1[-1]  # index `0～num-1` is the SWT data of different series
        >>> save_data_arr2 = data_arr2[-1]  # index `0～num-1` is the SWT data of different series
        >>> # You can save SWT to audio file
        >>> # af.write('SAVE_PATH_1.wav', save_data_arr1)
        >>> # af.write('SAVE_PATH_2.wav', save_data_arr2)

        Show wave plot

        >>> import matplotlib.pyplot as plt
        >>> from audioflux.display import fill_wave
        >>>
        >>> n_num, n_time = data_arr2.shape
        >>> fig, ax = plt.subplots(nrows=n_num + 1, figsize=(8, 16), sharex=True, layout='tight')
        >>> _ax0 = fill_wave(audio_arr[:n_time], samplate=sr, axes=ax[0])
        >>> _ax0.set_title('SWT arr2 origin')
        >>> for i in range(n_num):
        >>>     _ax = fill_wave(data_arr2[i], samplate=sr, axes=ax[i + 1])
        >>>     _ax.set_title(f'SWT arr2 num={i}')
        >>>
        >>> n_num, n_time = data_arr1.shape
        >>> fig, ax = plt.subplots(nrows=n_num + 1, figsize=(8, 16), sharex=True, layout='tight')
        >>> _ax0 = fill_wave(audio_arr[:n_time], samplate=sr, axes=ax[0])
        >>> _ax0.set_title('SWT arr1 origin')
        >>> for i in range(n_num):
        >>>     _ax = fill_wave(data_arr1[i], samplate=sr, axes=ax[i + 1])
        >>>     _ax.set_title(f'SWT arr1 num={i}')
    """

    def __init__(self, num, fft_length, wavelet_type=WaveletDiscreteType.SYM, t1=4, t2=0):
        super(SWT, self).__init__(pointer(OpaqueSWT()))

        self.num = num
        self.fft_length = fft_length
        self.wavelet_type = wavelet_type
        self.t1 = t1
        self.t2 = t2

        fn = self._lib['swtObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueSWT)),
                       c_int, c_int,
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.fft_length),
           pointer(c_int(self.wavelet_type.value)),
           pointer(c_int(self.t1)),
           pointer(c_int(self.t2)))
        self._is_created = True

    def swt(self, data_arr):
        """
        Get swt matrix

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data

        Returns
        -------
        m_data_arr1: np.ndarray [shape=(..., fre, time)]
        m_data_arr2: np.ndarray [shape=(..., fre, time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['swtObj_swt']
        fn.argtypes = [POINTER(OpaqueSWT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        if data_arr.ndim == 1:
            m_data_arr1 = np.zeros((self.num, self.fft_length), dtype=np.float32)
            m_data_arr2 = np.zeros((self.num, self.fft_length), dtype=np.float32)
            fn(self._obj, data_arr, m_data_arr1, m_data_arr2)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            m_data_arr1 = np.zeros((channel_num, self.num, self.fft_length), dtype=np.float32)
            m_data_arr2 = np.zeros((channel_num, self.num, self.fft_length), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], m_data_arr1[i], m_data_arr2[i])
            m_data_arr1 = revoke_channel(m_data_arr1, o_channel_shape, 2)
            m_data_arr2 = revoke_channel(m_data_arr2, o_channel_shape, 2)
        return m_data_arr1, m_data_arr2

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['swtObj_free']
            free_fn.argtypes = [POINTER(OpaqueSWT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
