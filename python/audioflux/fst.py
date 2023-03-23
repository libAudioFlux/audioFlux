import warnings

import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.base import Base
from audioflux.utils import check_audio, check_audio_length, format_channel, revoke_channel

__all__ = ['FST']


class OpaqueFST(Structure):
    _fields_ = []


class FST(Base):
    """
    Fast S-Transform (ST)

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    min_index: int
        Min bank index.

        Available range: (0, max_index)

    max_index: int
        Max bank index.

        Available range: (min_index, fft_length/2]

    See Also
    --------
    CQT
    ST
    DWT
    WPT
    SWT

    Examples
    --------

    Read 880Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('880')
    >>> audio_arr, sr = af.read(audio_path)
    >>> # FST can only input fft_length data
    >>> # For radix2_exp=12, then fft_length=4096
    >>> audio_arr = audio_arr[..., :4096]

    Create FST object

    >>> min_index, max_index = 1, 1024    # frequency is about 7.8125~8000Hz
    >>> obj = af.FST(radix2_exp=12, min_index=min_index, max_index=max_index, samplate=sr)

    Extract spectrogram

    >>> import numpy as np
    >>> spec_arr = obj.fst(audio_arr)
    >>> spec_arr = np.abs(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='FST Spectrogram')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, radix2_exp=12, min_index=1, max_index=None, samplate=32000):
        super(FST, self).__init__(pointer(OpaqueFST()))

        fft_length = 1 << radix2_exp

        if max_index is None:
            max_index = (fft_length // 2) - 1

        if min_index < 1:
            raise ValueError(f'min_index={min_index} must be a positive integer.')
        if max_index >= (fft_length / 2):
            raise ValueError(f'max_index={max_index} must be less than or equal to fft_length/2={fft_length / 2}')
        if min_index >= max_index:
            raise ValueError(f'min_index={min_index} must be less than max_index={max_index}')

        self.radix2_exp = radix2_exp
        self.min_index = min_index
        self.max_index = max_index
        self.samplate = samplate

        self.fft_length = fft_length
        self.num = self.max_index - self.min_index + 1

        fn = self._lib['fstObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueFST)), c_int]
        fn(self._obj,
           c_int(self.radix2_exp))
        self._is_created = True

    def get_fre_band_arr(self):
        """
        Get an array of FST frequency bands of different scales.

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """

        return np.arange(self.min_index, self.max_index + 1, dtype=np.float32) * self.samplate / self.fft_length

    def fst(self, data_arr):
        """
        Get spectrogram data

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., 2**radix2_exp)]
            Input audio data

        Returns
        -------
        out: np.ndarray [shape=(..., fre, time), dtype=np.complex]
        """

        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)
        data_arr = check_audio_length(data_arr, self.radix2_exp)

        st_fn = self._lib['fstObj_fst']
        st_fn.argtypes = [POINTER(OpaqueFST),
                          np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                          c_int,
                          c_int,
                          np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                          np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                          ]

        if data_arr.ndim == 1:
            m_real_arr = np.zeros((self.num, self.fft_length), dtype=np.float32)
            m_imag_arr = np.zeros((self.num, self.fft_length), dtype=np.float32)
            st_fn(self._obj, data_arr, c_int(self.min_index), c_int(self.max_index), m_real_arr, m_imag_arr)
            m_ret_arr = m_real_arr + (m_imag_arr * 1j)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            m_real_arr = np.zeros((channel_num, self.num, self.fft_length), dtype=np.float32)
            m_imag_arr = np.zeros((channel_num, self.num, self.fft_length), dtype=np.float32)
            for i in range(channel_num):
                st_fn(self._obj, data_arr[i], c_int(self.min_index), c_int(self.max_index),
                      m_real_arr[i], m_imag_arr[i])
            m_ret_arr = m_real_arr + (m_imag_arr * 1j)
            m_ret_arr = revoke_channel(m_ret_arr, o_channel_shape, 2)

        return m_ret_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        fre_band_arr = self.get_fre_band_arr()
        y_coords = np.insert(fre_band_arr, 0, fre_band_arr[0])
        return y_coords

    def x_coords(self):
        """
        Get the X-axis coordinate

        Returns
        -------
        out: [shape=(time,)]
        """
        x_coords = np.linspace(0, self.fft_length / self.samplate, self.fft_length + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['fstObj_free']
            free_fn.argtypes = [POINTER(OpaqueFST)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
