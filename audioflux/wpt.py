import warnings
import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p
from audioflux.type import WaveletDiscreteType
from audioflux.base import Base
from audioflux.utils import check_audio, check_audio_length

__all__ = ["WPT"]


class OpaqueWPT(Structure):
    _fields_ = []


class WPT(Base):
    """
    Wave Packet Transform (WPT)

    Parameters
    ----------
    num: int or None
        Number of frequency bins is ``2**num``

        If num is `None`, then ``num = radix2_exp - 1``

    radix2_exp: int
        ``fft_length=2**radix2_exp``

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
    SWT

    Examples
    --------

    Read 880Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('880')
    >>> audio_arr, sr = af.read(audio_path)
    >>> audio_arr = audio_arr[:4096]
    array([-5.5879354e-09, -9.3132257e-09,  0.0000000e+00, ...,
           -1.3013066e-01, -8.8289484e-02, -7.2360471e-02], dtype=float32)

    Create WPT object

    >>> from audioflux.type import WaveletDiscreteType
    >>> obj = af.WPT(num=7, radix2_exp=12, samplate=sr,
    >>>              wavelet_type=WaveletDiscreteType.SYM,
    >>>              t1=4, t2=0)

    Extract WPT data

    >>> import numpy as np
    >>> coef_arr, m_data_arr = obj.wpt(audio_arr)
    >>> m_data_arr = np.abs(m_data_arr)
    >>> coef_arr
    array([ 3.6137632e-01, -2.3430534e-01, -2.9199639e-01, ...,
            4.3719047e-06,  5.2029663e-04, -7.2287242e-03], dtype=float32)

    >>> m_data_arr
    array([[0.36137632, 0.36137632, 0.36137632, ..., 0.1472302 , 0.1472302 ,
            0.1472302 ],
           [0.07282908, 0.07282908, 0.07282908, ..., 0.07150738, 0.07150738,
            0.07150738],
           [0.03349842, 0.03349842, 0.03349842, ..., 0.01227883, 0.01227883,
            0.01227883],
           ...,
           [0.00414347, 0.00414347, 0.00414347, ..., 0.00337636, 0.00337636,
            0.00337636],
           [0.00756872, 0.00756872, 0.00756872, ..., 0.00425741, 0.00425741,
            0.00425741],
           [0.00048361, 0.00048361, 0.00048361, ..., 0.00722872, 0.00722872,
            0.00722872]], dtype=float32)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec, fill_plot
    >>> fig, ax = plt.subplots(nrows=2, sharex=True)
    >>> fill_spec(m_data_arr, axes=ax[0],
    >>>           x_coords=obj.x_coords(), y_coords=obj.y_coords(),
    >>>           x_axis='time', y_axis='log',
    >>>           title='WPT-abs')
    >>> fill_plot(obj.x_coords()[1:], coef_arr,
    >>>           axes=ax[1], label='WPT-coef')

    """

    def __init__(self, num=None, radix2_exp=12, samplate=32000,
                 wavelet_type=WaveletDiscreteType.SYM, t1=4, t2=0):
        super(WPT, self).__init__(pointer(OpaqueWPT()))

        if num is None:
            num = radix2_exp - 1

        # if num >= radix2_exp or num <= 0:
        #     raise ValueError(f'The num={num} range is [1, {radix2_exp - 1}]')

        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.wavelet_type = wavelet_type
        self.t1 = t1
        self.t2 = t2

        self.fft_length = 1 << self.radix2_exp

        fn = self._lib['wptObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueWPT)),
                       c_int, c_int,
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           c_int(self.num),
           c_int(self.radix2_exp),
           pointer(c_int(self.wavelet_type.value)),
           pointer(c_int(self.t1)),
           pointer(c_int(self.t2)))
        self._is_created = True

    def get_fre_band_arr(self):
        return np.linspace(0, 16000, (1 << self.num), dtype=np.float32)

    def wpt(self, data_arr):
        """
        Get wpt matrix

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)]
            Input audio data

        Returns
        -------
        out: np.ndarray [shape=(fre, time)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr)
        data_arr = check_audio_length(data_arr, self.radix2_exp)

        fn = self._lib['wptObj_wpt']
        fn.argtypes = [POINTER(OpaqueWPT),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                       ]

        coef_arr = np.zeros(self.fft_length, dtype=np.float32)
        m_data_arr = np.zeros((1 << self.num, self.fft_length), dtype=np.float32)

        fn(self._obj, data_arr, coef_arr, m_data_arr)
        return coef_arr, m_data_arr

    def y_coords(self):
        """
        Get the Y-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(fre,)]
        """
        fre_arr = self.get_fre_band_arr()
        y_coords = np.insert(fre_arr, 0, fre_arr[0])
        return y_coords

    def x_coords(self):
        """
        Get the X-axis coordinate

        Returns
        -------
        out: np.ndarray [shape=(time,)]
        """
        x_coords = np.linspace(0, self.fft_length / self.samplate, self.fft_length + 1)
        return x_coords

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['wptObj_free']
            free_fn.argtypes = [POINTER(OpaqueWPT)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
