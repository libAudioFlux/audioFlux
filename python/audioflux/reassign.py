import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.type import ReassignType, WindowType
from audioflux.base import Base
from audioflux.utils import check_audio, ascontiguous_swapaxex, format_channel, revoke_channel

__all__ = ['Reassign']


class OpaqueReassign(Structure):
    _fields_ = []


class Reassign(Base):
    """
    Reassign - reassign transform for STFT

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

    re_type: ReassignType
        Reassign type

    thresh: float
        thresh

    is_padding: bool
        Whether to use padding

    See Also
    --------
    Synsq
    WSST

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Create Reassign object

    >>> from audioflux.type import ReassignType, WindowType
    >>> obj = af.Reassign(radix2_exp=12, samplate=sr, window_type=WindowType.HANN,
    >>>                   slide_length=None, re_type=ReassignType.ALL)

    Extract spectrogram

    >>> import numpy as np
    >>> re_spec_arr, bft_spec_arr = obj.reassign(audio_arr)
    >>> re_spec_arr = np.abs(re_spec_arr)
    >>> bft_spec_arr = np.abs(bft_spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> audio_len = audio_arr.shape[-1]
    >>> # Show BFT/STFT Spectrogram
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(bft_spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='BFT/STFT Spectrogram')
    >>> fig.colorbar(img, ax=ax)
    >>>
    >>> # Show Reassign Spectrogram
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(re_spec_arr, axes=ax,
    >>>                 x_coords=obj.x_coords(audio_len),
    >>>                 y_coords=obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Reassign Spectrogram')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, radix2_exp=12, samplate=32000, window_type=WindowType.HANN,
                 slide_length=None, re_type=ReassignType.ALL, thresh=0.001,
                 is_padding=False):

        super(Reassign, self).__init__(pointer(OpaqueReassign()))

        self.fft_length = fft_length = 1 << radix2_exp

        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.window_type = window_type
        if slide_length is None:
            slide_length = fft_length // 4
        self.slide_length = slide_length
        self.re_type = re_type
        self.thresh = thresh
        self.is_padding = is_padding
        self.is_continue = False
        self.order = 1
        self.result_type = 0

        fn = self._lib['reassignObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueReassign)),
                       c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int),
                       POINTER(c_int), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int)]
        fn(self._obj,
           c_int(self.radix2_exp),
           pointer(c_int(self.samplate)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.re_type.value)),
           pointer(c_float(self.thresh)),
           pointer(c_int(int(self.is_padding))),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True

    def cal_time_length(self, data_length):
        """
        Calculate the length of a frame from audio data.

        Parameters
        ----------
        data_length: int
            The length of the data to be calculated.

        Returns
        -------
        out: int
        """
        c_fn = self._lib['reassignObj_calTimeLength']
        c_fn.argtypes = [POINTER(OpaqueReassign), c_int]
        return c_fn(self._obj, c_int(data_length))

    def set_result_type(self, result_type):
        """
        Set result type.

        Parameters
        ----------
        result_type: int, 0 or 1
            - If `0`, then the result is a matrix of complex numbers.
            - If `1`, then the result is a matrix of real numbers.
        """
        c_fn = self._lib['reassignObj_setResultType']
        c_fn.argtypes = [POINTER(OpaqueReassign), c_int]
        c_fn(self._obj, c_int(result_type))
        self.result_type = result_type

    def set_order(self, order):
        """
        Set order

        Parameters
        ----------
        order: int
            order >= 1
        """
        c_fn = self._lib['reassignObj_setOrder']
        c_fn.argtypes = [POINTER(OpaqueReassign), c_int]
        c_fn(self._obj, c_int(order))
        self.order = order

    def reassign(self, data_arr, result_type=0):
        """
        Get reassign matrix

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio data

        result_type: int, 0 or 1
            - If `0`, then the result is a matrix of complex numbers.
            - If `1`, then the result is a matrix of real numbers.

        Returns
        -------
        m_arr1: np.ndarray [shape=(..., fre, time), dtype=(np.complex or np.float32)]
            The matrix of reassign

        m_arr2: np.ndarray [shape=(..., fre, time), dtype=np.complex]
            The matrix of origin(BFT/STFT)
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        if result_type != self.result_type:
            self.set_result_type(result_type)

        c_fn = self._lib['reassignObj_reassign']
        c_fn.argtypes = [POINTER(OpaqueReassign),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                         c_int,
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')]

        data_len = data_arr.shape[-1]
        n_time = self.cal_time_length(data_len)
        n_fre = self.fft_length // 2 + 1

        if data_arr.ndim == 1:
            shape = (n_time, n_fre)
            m_real_arr1 = np.zeros(shape, dtype=np.float32)
            m_imag_arr1 = np.zeros(shape, dtype=np.float32)
            m_real_arr2 = np.zeros(shape, dtype=np.float32)
            m_imag_arr2 = np.zeros(shape, dtype=np.float32)
            c_fn(self._obj, data_arr, c_int(data_len), m_real_arr1, m_imag_arr1, m_real_arr2, m_imag_arr2)
            m_arr1 = (m_real_arr1 + m_imag_arr1 * 1j) if self.result_type == 0 else m_real_arr1
            m_arr2 = m_real_arr2 + m_imag_arr2 * 1j
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            shape = (channel_num, n_time, n_fre)
            m_real_arr1 = np.zeros(shape, dtype=np.float32)
            m_imag_arr1 = np.zeros(shape, dtype=np.float32)
            m_real_arr2 = np.zeros(shape, dtype=np.float32)
            m_imag_arr2 = np.zeros(shape, dtype=np.float32)

            for i in range(channel_num):
                c_fn(self._obj, data_arr[i], c_int(data_len),
                     m_real_arr1[i], m_imag_arr1[i], m_real_arr2[i], m_imag_arr2[i])
            m_arr1 = (m_real_arr1 + m_imag_arr1 * 1j) if self.result_type == 0 else m_real_arr1
            m_arr2 = m_real_arr2 + m_imag_arr2 * 1j
            m_arr1 = revoke_channel(m_arr1, o_channel_shape, 2)
            m_arr2 = revoke_channel(m_arr2, o_channel_shape, 2)

        m_arr1 = ascontiguous_swapaxex(m_arr1, -1, -2)
        m_arr2 = ascontiguous_swapaxex(m_arr2, -1, -2)
        return m_arr1, m_arr2

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        m_len = self.fft_length // 2 + 1 + 1
        y_coords = np.linspace(0, self.samplate // 2, m_len)
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
        x_coords = np.linspace(0, data_length * 1. / self.samplate,
                               self.cal_time_length(data_length) + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['reassignObj_free']
            free_fn.argtypes = [POINTER(OpaqueReassign)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
