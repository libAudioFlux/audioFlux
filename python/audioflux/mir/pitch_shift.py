import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel

__all__ = ["PitchShift"]


class OpaquePitchShift(Structure):
    _fields_ = []


class PitchShift(Base):
    """
    Pitch shift algorithm

    Parameters
    ----------
    radix2_exp: int
        ``fft_length=2**radix2_exp``

    slide_length: int
        Window sliding length.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Examples
    --------

    Read voice audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('voice')
    >>> audio_arr, sr = af.read(audio_path)

    Compute pitch shift

    >>> pitch_shift_obj = af.PitchShift(radix2_exp=12, slide_length=1024, window_type=af.type.WindowType.HANN)
    >>> audio_arr1 = pitch_shift_obj.pitch_shift(audio_arr, n_semitone=6, samplate=sr)
    >>> # af.write('./audio_arr1.wav', audio_arr1, sr)

    Show plot

    >>> yin_obj = af.PitchYIN(samplate=sr)
    >>> pitch_arr, _, _ = yin_obj.pitch(audio_arr)
    >>> pitch_arr1, _, _ = yin_obj.pitch(audio_arr1)

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> times = np.arange(len(pitch_arr)) * (yin_obj.slide_length / sr)
    >>> ax.scatter(times, pitch_arr, s=2, label='Original')
    >>> ax.scatter(times, pitch_arr1, s=2, label='Pitch Shift')
    >>> ax.set_ylim(200, 750)
    >>> ax.legend()
    """
    def __init__(self, radix2_exp=12, slide_length=1024, window_type=WindowType.HANN):
        super(PitchShift, self).__init__(pointer(OpaquePitchShift()))

        self.radix2_exp = radix2_exp
        self.window_type = window_type
        self.slide_length = slide_length

        fn = self._lib['pitchShiftObj_new']
        fn.argtypes = [POINTER(POINTER(OpaquePitchShift)),
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]

        fn(self._obj,
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.window_type.value)))
        self._is_created = True

    def pitch_shift(self, data_arr, n_semitone, samplate=32000):
        """
        Compute the time stretch

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Audio data array

        n_semitone: int
            Pitch shift in semitone. -12 <= n_semitone <= 12.

        samplate: int
            Sample rate

        Returns
        -------
        arr: np.ndarray [shape=(..., n)]
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        if n_semitone < -12 or n_semitone > 12:
            raise ValueError(f"n_semitone={n_semitone} must be in range [-12, 12]")

        fn = self._lib['pitchShiftObj_pitchShift']
        fn.argtypes = [
            POINTER(OpaquePitchShift),
            c_int,
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int,
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        ]

        data_len = data_arr.shape[-1]
        data_len_c = c_int(data_len)
        n_semitone_c = c_int(n_semitone)
        samplate_c = c_int(samplate)

        if data_arr.ndim == 1:
            new_arr = np.zeros(data_len, dtype=np.float32)
            fn(self._obj, samplate_c, n_semitone_c, data_arr, data_len_c, new_arr)
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            new_arr = np.zeros((channel_num, data_len), dtype=np.float32)
            for i in range(channel_num):
                fn(self._obj, samplate_c, n_semitone_c, data_arr[i], data_len_c, new_arr[i])

            new_arr = revoke_channel(new_arr, o_channel_shape, 1)

        return new_arr

    # def __del__(self):
    #     if self._is_created:
    #         fn = self._lib['pitchShiftObj__free']
    #         fn.argtypes = [POINTER(OpaquePitchShift)]
    #         fn.restype = c_void_p
    #         fn(self._obj)
