import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_float, c_void_p
from audioflux.base import Base
from audioflux.utils import format_channel, revoke_channel

__all__ = ['Synsq']


class OpaqueSynsq(Structure):
    _fields_ = []


class Synsq(Base):
    """
    Synsq algorithm

    Parameters
    ----------
    num: int
        Number of frequency bins to generate

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio

    order: int
        order value

    thresh: float
        thresh value

    See Also
    --------
    Reassign
    WSST

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)
    >>> # WSST can only input fft_length data
    >>> # For radix2_exp=12, then fft_length=4096
    >>> audio_arr = audio_arr[..., :4096]

    Create CWT object of octave

    >>> from audioflux.type import SpectralFilterBankScaleType, WaveletContinueType
    >>> from audioflux.utils import note_to_hz

    >>> cwt_obj = af.CWT(num=84, radix2_exp=12, samplate=sr, low_fre=note_to_hz('C1'),
    >>>                  bin_per_octave=12, wavelet_type=WaveletContinueType.MORLET,
    >>>                  scale_type=SpectralFilterBankScaleType.OCTAVE)

    Extract CWT spectrogram

    >>> import numpy as np
    >>> cwt_spec_arr = cwt_obj.cwt(audio_arr)

    Create Synsq object

    >>> synsq_obj = af.Synsq(num=cwt_obj.num,
    >>>                      radix2_exp=cwt_obj.radix2_exp,
    >>>                      samplate=cwt_obj.samplate)

    Extract Synsq data

    >>> synsq_arr = synsq_obj.synsq(cwt_spec_arr,
    >>>                             filter_bank_type=cwt_obj.scale_type,
    >>>                             fre_arr=cwt_obj.get_fre_band_arr())

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> # Show CWT
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(np.abs(cwt_spec_arr), axes=ax,
    >>>                 x_coords=cwt_obj.x_coords(),
    >>>                 y_coords=cwt_obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='CWT')
    >>> fig.colorbar(img, ax=ax)
    >>>
    >>> # Show Synsq
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(np.abs(synsq_arr), axes=ax,
    >>>                 x_coords=cwt_obj.x_coords(),
    >>>                 y_coords=cwt_obj.y_coords(),
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Synsq')
    >>> fig.colorbar(img, ax=ax)
    """

    def __init__(self, num, radix2_exp=12, samplate=32000, order=1, thresh=0.001):
        super(Synsq, self).__init__(pointer(OpaqueSynsq()))

        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.order = order
        self.thresh = thresh

        fn = self._lib['synsqObj_new']
        fn.argtypes = [
            POINTER(POINTER(OpaqueSynsq)),
            c_int,
            c_int,
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_float),
        ]
        fn(self._obj,
           c_int(self.num),
           c_int(self.radix2_exp),
           pointer(c_int(self.samplate)),
           pointer(c_int(self.order)),
           pointer(c_float(self.thresh)))
        self._is_created = True

    def synsq(self, m_data_arr, filter_bank_type, fre_arr):
        """
        Get synsq data

        Parameters
        ----------
        m_data_arr: np.ndarray [shape=(..., fre, time), dtype=np.complex]
            spectrogram data

        filter_bank_type: SpectralFilterBankScaleType
            bank type，need to correspond to `m_data_arr`.

            See: `type.SpectralFilterBankScaleType`

        fre_arr: np.ndarray [shape=(fre,)]
            fre band arr，need to correspond to `m_data_arr`.
            Use `get_fre_band_arr` method to get

        Returns
        -------
        out: np.ndarray [shape=(..., fre, time), dtype=np.complex]
        """

        if not np.iscomplexobj(m_data_arr):
            raise ValueError(f"m_data_arr with dtype={m_data_arr.dtype} is not of complex type")
        if m_data_arr.ndim < 2:
            raise ValueError(f"m_data_arr.ndim=[{m_data_arr.ndim}] should be greater than 1")

        c_fn = self._lib['synsqObj_synsq']
        c_fn.argtypes = [POINTER(OpaqueSynsq),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                         c_int,
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
                         np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS')]

        if m_data_arr.ndim == 2:
            m_real_arr1 = np.ascontiguousarray(m_data_arr.real)
            m_imag_arr1 = np.ascontiguousarray(m_data_arr.imag)
            m_real_arr2 = np.zeros_like(m_data_arr, dtype=np.float32)
            m_imag_arr2 = np.zeros_like(m_data_arr, dtype=np.float32)
            c_fn(self._obj, fre_arr, c_int(filter_bank_type.value),
                 m_real_arr1, m_imag_arr1, m_real_arr2, m_imag_arr2)
            m_synsq_arr = m_real_arr2 + m_imag_arr2 * 1j
        else:
            m_data_arr, o_channel_shape = format_channel(m_data_arr, 2)
            channel_num = m_data_arr.shape[0]

            m_real_arr1 = np.ascontiguousarray(m_data_arr.real)
            m_imag_arr1 = np.ascontiguousarray(m_data_arr.imag)
            m_real_arr2 = np.zeros_like(m_real_arr1, dtype=np.float32)
            m_imag_arr2 = np.zeros_like(m_imag_arr1, dtype=np.float32)
            for i in range(channel_num):
                c_fn(self._obj, fre_arr, c_int(filter_bank_type.value),
                     m_real_arr1[i], m_imag_arr1[i], m_real_arr2[i], m_imag_arr2[i])
            m_synsq_arr = m_real_arr2 + m_imag_arr2 * 1j
            m_synsq_arr = revoke_channel(m_synsq_arr, o_channel_shape, 2)

        return m_synsq_arr

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['synsqObj_free']
            free_fn.argtypes = [POINTER(OpaqueSynsq)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
