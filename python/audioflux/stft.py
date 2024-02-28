import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.type import WindowType, PaddingPositionType, PaddingModeType
from audioflux.base import Base
from audioflux.utils import check_audio, format_channel, revoke_channel, ascontiguous_swapaxex

__all__ = ["STFT", 'OpaqueSTFT']


class OpaqueSTFT(Structure):
    _fields_ = []


class STFT(Base):
    """
    Short-time Fourier transform (STFT).

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    See Also
    --------
    BFT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Compute stft and istft

    >>> stft_obj = af.STFT(radix2_exp=12, window_type=af.type.WindowType.RECT, slide_length=1024)
    >>> spec_arr = stft_obj.stft(audio_arr)
    >>> new_audio_arr = stft_obj.istft(spec_arr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_title('STFT Spectrogram')
    >>> img = af.display.fill_spec(np.abs(spec_arr), axes=ax,
    >>>                            y_coords=stft_obj.y_coords(sr),
    >>>                            x_coords=stft_obj.x_coords(audio_arr.shape[-1], sr),
    >>>                            y_axis='log', x_axis='time')
    >>> fig.colorbar(img, ax=ax)

    >>> fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)
    >>> ax = af.display.fill_wave(audio_arr, axes=axes[0])
    >>> ax.set_title('Original')
    >>> ax = af.display.fill_wave(new_audio_arr, axes=axes[1])
    >>> ax.set_title('ISTFT Result')

    """

    def __init__(self, radix2_exp=12, window_type=WindowType.RECT, slide_length=1024):
        super(STFT, self).__init__(pointer(OpaqueSTFT()))

        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length
        self.is_continue = False

        self.is_pad = False
        self.position_type = PaddingPositionType.CENTER
        self.mode_type = PaddingModeType.CONSTANT
        self.pad_value1 = 0.
        self.pad_value2 = 0.

        self.fft_length = 1 << radix2_exp

        fn = self._lib['stftObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueSTFT)),
                       c_int, POINTER(c_int),
                       POINTER(c_int), POINTER(c_int)]
        fn(self._obj,
           c_int(self.radix2_exp),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True

    def set_slide_length(self, slide_length):
        """
        Set the slide length.

        Parameters
        ----------
        slide_length: int
            Window sliding length.
        """
        fn = self._lib['stftObj_setSlideLength']
        fn.argtypes = [POINTER(OpaqueSTFT), c_int]
        fn(self._obj, c_int(slide_length))
        self.slide_length = slide_length

    def set_padding(self, position_type=PaddingPositionType.CENTER, mode_type=PaddingModeType.CONSTANT,
                    value1=0.0, value2=0.0):
        """
        Set padding parameters.

        Before calling this function, you must use enable_padding and set is_pad to True.

        Parameters
        ----------
        position_type: PaddingPositionType
            Padding position type.

            See: `type.PaddingPositionType`

        mode_type: PaddingModeType
            Padding mode type.

            See: `type.PaddingModeType`

        value1: float
            Padding value1.

            If mode_type is `CONSTANT` and position_type is `CENTER`, value1 is the left padding value, length is `fft_length // 2`.
            If mode_type is `CONSTANT` and position_type is `LEFT`, value1 is the left padding value, length is `fft_length // 2`.
            If mode_type is `CONSTANT` and position_type is `RIGHT`, value1 is the right padding value, length is `fft_length // 2`.
            Other modes are not used.

        value2: float
            Padding value2.

            If mode_type is `CONSTANT` and position_type is `CENTER`, value2 is the right padding value, length is `fft_length // 2`.
            Other modes are not used.

        Returns
        -------

        """
        fn = self._lib['stftObj_setPadding']
        fn.argtypes = [POINTER(OpaqueSTFT),
                       POINTER(c_int), POINTER(c_int),
                       POINTER(c_float), POINTER(c_float)]
        fn(self._obj, pointer(c_int(position_type.value)), pointer(c_int(mode_type.value)),
           pointer(c_float(value1)), pointer(c_float(value2)))
        self.position_type = position_type
        self.mode_type = mode_type
        self.pad_value1 = value1
        self.pad_value2 = value2

    def use_window_data_arr(self, data_arr):
        """
        Custom window data array.

        Default window data array is generated by `window_type`.

        Parameters
        ----------
        data_arr: np.ndarray [shape=(fft_length,)]
            Window data array.

        Returns
        -------

        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        if data_arr.ndim != 1:
            raise ValueError(f'data_arr.ndim must be 1')
        if data_arr.shape[-1] != self.fft_length:
            raise ValueError(f'data_arr length[{data_arr.shape[-1]}] must be {self.fft_length}')

        fn = self._lib['stftObj_useWindowDataArr']
        fn.argtypes = [POINTER(OpaqueSTFT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')]
        fn(self._obj, data_arr)

    def get_window_data_arr(self):
        """
        Get window data array.

        Returns
        -------
        out: np.ndarray [shape=(fft_length,)]
            Window data array.
        """
        fn = self._lib['stftObj_getWindowDataArr']
        fn.argtypes = [POINTER(OpaqueSTFT)]
        fn.restype = c_void_p
        p = fn(self._obj)
        ret = np.frombuffer((c_float * self.fft_length).from_address(p), np.float32).copy()
        return ret

    def enable_padding(self, flag=False):
        """
        Whether to enable padding.

        Default is `False`.

        Parameters
        ----------
        flag: bool
            Whether to enable padding.

        Returns
        -------

        """
        fn = self._lib['stftObj_enablePadding']
        fn.argtypes = [POINTER(OpaqueSTFT), c_int]
        fn(self._obj, c_int(int(flag)))
        self.is_pad = flag

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        - ``fft_length = 2 ** radix2_exp``
        - ``(data_length - fft_length) // slide_length + 1``

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['stftObj_calTimeLength']
        fn.argtypes = [POINTER(OpaqueSTFT), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def cal_data_length(self, time_length):
        """
        Calculate the length of the audio data from the frame length.

        Parameters
        ----------
        time_length: int
            The length of the frame to be calculated.

        Returns
        -------
        out: int
        """
        fn = self._lib['stftObj_calDataLength']
        fn.argtypes = [POINTER(OpaqueSTFT), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(time_length))

    def stft(self, data_arr):
        """
        Calculate STFT (Short-Time Fourier Transform) data.

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data

        Returns
        -------
        out: np.ndarray [shape=(..., fft_length // 2 + 1, time_length), dtype=np.complex64]
            STFT data
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['stftObj_stft']
        fn.argtypes = [POINTER(OpaqueSTFT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        data_length = data_arr.shape[-1]
        time_length = self.cal_time_length(data_length)

        if data_arr.ndim == 1:
            m_real_arr = np.zeros((time_length, self.fft_length), dtype=np.float32)
            m_imag_arr = np.zeros((time_length, self.fft_length), dtype=np.float32)
            fn(self._obj, data_arr, data_length, m_real_arr, m_imag_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            m_real_arr = np.zeros((channel_num, time_length, self.fft_length), dtype=np.float32)
            m_imag_arr = np.zeros((channel_num, time_length, self.fft_length), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, data_arr[i], data_length, m_real_arr[i], m_imag_arr[i])
            m_real_arr = revoke_channel(m_real_arr, o_channel_shape, 2)
            m_imag_arr = revoke_channel(m_imag_arr, o_channel_shape, 2)
        m_data_arr = (m_real_arr + m_imag_arr * 1j)
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_data_arr = m_data_arr[..., :(self.fft_length // 2 + 1), :]
        return m_data_arr

    def istft(self, m_data_arr, method_type=0):
        """
        Calculate ISTFT (Inverse Short-Time Fourier Transform) data.

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fft_length // 2 + 1, time_length), dtype=np.complex64]
            Input STFT data

        method_type: int
            0: weight(default)
            1: overlap-add

        Returns
        -------
        out: np.ndarray [shape=(..., data_length), dtype=np.float32]
            ISTFT data
        """
        if not np.iscomplexobj(m_data_arr):
            raise ValueError(f'm_data_arr must be of type np.complex')
        if m_data_arr.ndim < 2:
            raise ValueError(f"m_data_arr's dimensions must be greater than 1")

        fn = self._lib['stftObj_istft']
        fn.argtypes = [POINTER(OpaqueSTFT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       c_int,
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        # fill fft_length // 2 + 1 => fft_length
        stft_arr1 = np.copy(m_data_arr)[..., ::-1, :][..., 1:-1, :]
        stft_arr1.imag *= -1
        m_data_arr = np.concatenate([m_data_arr, stft_arr1], axis=-2)  # shape=(..., fft_length, time_length)

        time_length = m_data_arr.shape[-1]
        m_data_arr = ascontiguous_swapaxex(m_data_arr, -1, -2)
        m_real_arr = m_data_arr.real.astype(dtype=np.float32)
        m_imag_arr = m_data_arr.imag.astype(dtype=np.float32)

        data_length = self.cal_data_length(time_length)

        if m_data_arr.ndim == 2:
            ret_arr = np.zeros(data_length, dtype=np.float32)
            fn(self._obj, m_real_arr, m_imag_arr, c_int(time_length), c_int(method_type), ret_arr)
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num, t_len, f_len = m_data_arr.shape

            ret_arr = np.zeros((channel_num, data_length), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, m_real_arr[i], m_imag_arr[i], c_int(time_length), c_int(method_type), ret_arr[i])
            ret_arr = revoke_channel(ret_arr, o_channel_shape, 1)
        return ret_arr

    def y_coords(self, samplate=32000):
        """
        Get the Y-axis coordinate

        Parameters
        ----------
        samplate: int
            Sampling rate of the incoming audio.

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        y_coords = np.linspace(0, samplate // 2, self.fft_length // 2 + 1)
        y_coords = np.insert(y_coords, 0, 0)
        return y_coords

    def x_coords(self, data_length, samplate=32000):
        """
        Get the X-axis coordinate

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        samplate: int
            Sampling rate of the incoming audio.

        Returns
        -------
        out: np.ndarray [shape=(time,)]
        """
        if data_length < self.fft_length:
            raise ValueError(f'radix2_exp={self.radix2_exp}(fft_length={self.fft_length}) '
                             f'is too large for data_length={data_length}')
        x_coords = np.linspace(0, data_length / samplate, self.cal_time_length(data_length) + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['stftObj_free']
            free_fn.argtypes = [POINTER(OpaqueSTFT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
