import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_float
from audioflux.base import Base
from audioflux.type import WindowType
from audioflux.utils import check_audio, format_channel, revoke_channel, ascontiguous_T

__all__ = ["PitchFFP"]


class OpaquePitchFFP(Structure):
    _fields_ = [
        # ('isContinue', c_int),
        #
        # ('stftObj', POINTER(c_void_p)),
        #
        # ('fftLength', c_int),
        # ('slideLength', c_int),
        #
        # ('cutLength', c_int),
        # ('peakLength', c_int),
        #
        # ('minIndex', c_int),
        # ('maxIndex', c_int),
        #
        # ('timeLength', c_int),
    ]


class PitchFFP(Base):
    """
    Pitch FFP algorithm

    Parameters
    ----------
    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. Default is `32.0`.

    high_fre: float
        Highest frequency. Default is `2000.0`.

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

    Extract pitch

    >>> pitch_obj = af.PitchFFP(samplate=sr, radix2_exp=12, low_fre=27, high_fre=2000,
    >>>                         slide_length=1024, window_type=af.type.WindowType.HAMM)
    >>> fre_arr, db_arr, extra_data_dic = pitch_obj.pitch(
    >>>     audio_arr,
    >>>     has_corr_data=True,
    >>>     has_cut_data=True,
    >>>     has_flag_data=True,
    >>>     has_light_data=True,
    >>>     has_temporal_data=True)
    >>> cut_corr_arr, cut_db_arr, cut_height_arr, cut_length_arr = extra_data_dic['cut_data']
    >>> flag_arr, = extra_data_dic['flag_data']
    >>> light_arr, = extra_data_dic['light_data']
    >>> avg_temp_arr, max_temp_arr, percent_temp_arr = extra_data_dic['temporal_data']

    Show pitch plot

    >>> import matplotlib.pyplot as plt
    >>> times = np.arange(fre_arr.shape[-1]) * (pitch_obj.slide_length / sr)
    >>>
    >>> # Show pitch plot
    >>> fig, ax = plt.subplots(nrows=3, sharex=True)
    >>> ax[0].set_title('PitchFFP')
    >>> af.display.fill_wave(audio_arr, axes=ax[0])
    >>> ax[1].scatter(times, fre_arr, s=2)
    >>> ax[2].scatter(times, db_arr, s=2)
    >>>
    >>> # Show cut data plot
    >>> fig, ax = plt.subplots(nrows=4, sharex=True)
    >>> ax[0].set_title('FFP Corr Data')
    >>> af.display.fill_wave(audio_arr, axes=ax[0])
    >>> for i in range(cut_corr_arr.shape[0]):
    >>>     ax[1].scatter(times, cut_corr_arr[i], s=2)
    >>> for i in range(cut_db_arr.shape[0]):
    >>>     ax[2].scatter(times, cut_db_arr[i], s=2)
    >>> for i in range(cut_height_arr.shape[0]):
    >>>     ax[3].scatter(times, cut_height_arr[i], s=2)
    >>>
    >>> # Show flag and light data plot
    >>> fig, ax = plt.subplots(nrows=3, sharex=True)
    >>> ax[0].set_title('FFP Flag and Light Data')
    >>> af.display.fill_wave(audio_arr, axes=ax[0])
    >>> ax[1].scatter(times, flag_arr, s=2, label='Flag')
    >>> ax[1].legend()
    >>> ax[2].scatter(times, light_arr, s=2, label='Light')
    >>> ax[2].legend()
    >>>
    >>> # Show temporal data plot
    >>> fig, ax = plt.subplots(nrows=4, sharex=True)
    >>> ax[0].set_title('FFP Temporal Data')
    >>> af.display.fill_wave(audio_arr, axes=ax[0])
    >>> ax[1].scatter(times, avg_temp_arr, s=2)
    >>> ax[2].scatter(times, max_temp_arr, s=2)
    >>> ax[3].scatter(times, percent_temp_arr, s=2)
    """

    def __init__(self, samplate=32000, low_fre=32.0, high_fre=2000.0, radix2_exp=12,
                 slide_length=1024, window_type=WindowType.HAMM):
        super(PitchFFP, self).__init__(pointer(OpaquePitchFFP()))

        if low_fre >= high_fre:
            raise ValueError(f'`low_fre` must be smaller than `high_fre`')

        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.radix2_exp = radix2_exp
        self.slide_length = slide_length
        self.window_type = window_type
        self.is_continue = False

        self.temp_base = -18.0

        fn = self._lib['pitchFFPObj_new']
        fn.argtypes = [POINTER(POINTER(OpaquePitchFFP)),
                       POINTER(c_int), POINTER(c_float), POINTER(c_float),
                       POINTER(c_int), POINTER(c_int), POINTER(c_int),
                       POINTER(c_int)]

        fn(self._obj,
           pointer(c_int(self.samplate)),
           pointer(c_float(self.low_fre)),
           pointer(c_float(self.high_fre)),
           pointer(c_int(self.radix2_exp)),
           pointer(c_int(self.slide_length)),
           pointer(c_int(self.window_type.value)),
           pointer(c_int(int(self.is_continue))))
        self._is_created = True

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
        fn = self._lib['pitchFFPObj_calTimeLength']
        fn.argtypes = [POINTER(OpaquePitchFFP), c_int]
        fn.restype = c_int
        return fn(self._obj, c_int(data_length))

    def set_temp_base(self, temp_base):
        """
        Set temproal base

        Parameters
        ----------
        temp_base: float
            temproal base

        """
        if not -36 < temp_base < 0:
            raise ValueError(f'`temp_base` must be between -36 and 0.')

        fn = self._lib['pitchFFPObj_setTempBase']
        fn.argtypes = [POINTER(OpaquePitchFFP), c_float]
        fn(self._obj, c_float(temp_base))

        self.temp_base = temp_base

    def _pitch(self, data_arr):
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        fn = self._lib['pitchFFPObj_pitch']
        fn.argtypes = [POINTER(OpaquePitchFFP),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       c_int,
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
                       ]

        data_len = data_arr.shape[-1]
        time_length = self.cal_time_length(data_len)

        fre_arr = np.zeros(time_length, dtype=np.float32)
        db_arr = np.zeros(time_length, dtype=np.float32)
        fn(self._obj, data_arr, c_int(data_len), fre_arr, db_arr)

        return fre_arr, db_arr

    def _get_corr_data(self, data_length):
        fn = self._lib['pitchFFPObj_getCorrData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        pp_corr_arr = pointer(c_void_p())
        pp_db_arr = pointer(c_void_p())
        pp_height_arr = pointer(c_void_p())
        pp_len_arr = pointer(c_void_p())

        time_len = self.cal_time_length(data_length)

        print('xxxxx', time_len)

        peak_len = fn(self._obj, pp_corr_arr, pp_db_arr, pp_height_arr, pp_len_arr)

        corr_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_corr_arr[0]), np.float32).copy()
        db_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_db_arr[0]), np.float32).copy()
        height_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_height_arr[0]), np.float32).copy()
        len_arr = np.frombuffer((c_int * time_len).from_address(pp_len_arr[0]), np.int32).copy()

        corr_arr = ascontiguous_T(corr_arr.reshape((time_len, peak_len)))
        db_arr = ascontiguous_T(db_arr.reshape((time_len, peak_len)))
        height_arr = ascontiguous_T(height_arr.reshape((time_len, peak_len)))

        return corr_arr, db_arr, height_arr, len_arr

    def _get_cut_data(self, data_length):

        fn = self._lib['pitchFFPObj_getCutData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        pp_corr_arr = pointer(c_void_p())
        pp_db_arr = pointer(c_void_p())
        pp_height_arr = pointer(c_void_p())
        pp_len_arr = pointer(c_void_p())

        time_len = self.cal_time_length(data_length)

        peak_len = fn(self._obj, pp_corr_arr, pp_db_arr, pp_height_arr, pp_len_arr)

        corr_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_corr_arr[0]), np.float32).copy()
        db_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_db_arr[0]), np.float32).copy()
        height_arr = np.frombuffer((c_float * time_len * peak_len).from_address(pp_height_arr[0]), np.float32).copy()
        len_arr = np.frombuffer((c_int * time_len).from_address(pp_len_arr[0]), np.int32).copy()

        corr_arr = ascontiguous_T(corr_arr.reshape((time_len, peak_len))[:, :4])
        db_arr = ascontiguous_T(db_arr.reshape((time_len, peak_len))[:, :4])
        height_arr = ascontiguous_T(height_arr.reshape((time_len, peak_len))[:, :4])

        return corr_arr, db_arr, height_arr, len_arr

    def _get_flag_data(self):

        fn = self._lib['pitchFFPObj_getFlagData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        pp_flag_arr = pointer(c_void_p())
        time_len = fn(self._obj, pp_flag_arr)
        flag_arr = np.frombuffer((c_int * time_len).from_address(pp_flag_arr[0]), np.int32).copy()
        return flag_arr

    def _get_light_data(self):

        fn = self._lib['pitchFFPObj_getLightData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        pp_light_arr = pointer(c_void_p())
        time_len = fn(self._obj, pp_light_arr)
        light_arr = np.frombuffer((c_float * time_len).from_address(pp_light_arr[0]), np.float32).copy()
        return light_arr

    def _get_format_data(self, data_length):
        fn = self._lib['pitchFFPObj_getFormatData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        time_len = self.cal_time_length(data_length)

        pp_format_flag_arr = pointer(c_void_p())
        pp_fre_arr1 = pointer(c_void_p())
        pp_fre_arr2 = pointer(c_void_p())
        pp_fre_arr3 = pointer(c_void_p())
        pp_db_arr1 = pointer(c_void_p())
        pp_db_arr2 = pointer(c_void_p())
        pp_db_arr3 = pointer(c_void_p())

        fn(self._obj, pp_format_flag_arr,
           pp_fre_arr1, pp_fre_arr2, pp_fre_arr3,
           pp_db_arr1, pp_db_arr2, pp_db_arr3)

        format_flag_arr = np.frombuffer((c_int * time_len).from_address(pp_format_flag_arr[0]), np.int32).copy()
        fre_arr1 = np.frombuffer((c_float * time_len).from_address(pp_fre_arr1[0]), np.float32).copy()
        fre_arr2 = np.frombuffer((c_float * time_len).from_address(pp_fre_arr2[0]), np.float32).copy()
        fre_arr3 = np.frombuffer((c_float * time_len).from_address(pp_fre_arr3[0]), np.float32).copy()
        db_arr1 = np.frombuffer((c_float * time_len).from_address(pp_db_arr1[0]), np.float32).copy()
        db_arr2 = np.frombuffer((c_float * time_len).from_address(pp_db_arr2[0]), np.float32).copy()
        db_arr3 = np.frombuffer((c_float * time_len).from_address(pp_db_arr3[0]), np.float32).copy()

        return format_flag_arr, fre_arr1, fre_arr2, fre_arr3, db_arr1, db_arr2, db_arr3

    def _get_temporal_data(self):

        fn = self._lib['pitchFFPObj_getTemporalData']
        fn.argtypes = [
            POINTER(OpaquePitchFFP),
            POINTER(c_void_p),
            POINTER(c_void_p),
            POINTER(c_void_p),
        ]
        fn.restype = c_int

        pp_avg_temp_arr = pointer(c_void_p())
        pp_max_temp_arr = pointer(c_void_p())
        pp_percent_temp_arr = pointer(c_void_p())

        time_len = fn(self._obj, pp_avg_temp_arr, pp_max_temp_arr, pp_percent_temp_arr)

        avg_temp_arr = np.frombuffer((c_float * time_len).from_address(pp_avg_temp_arr[0]), np.float32).copy()
        max_temp_arr = np.frombuffer((c_float * time_len).from_address(pp_max_temp_arr[0]), np.float32).copy()
        percent_temp_arr = np.frombuffer((c_float * time_len).from_address(pp_percent_temp_arr[0]), np.float32).copy()

        return avg_temp_arr, max_temp_arr, percent_temp_arr

    def pitch(self, data_arr, has_corr_data=False, has_cut_data=False, has_flag_data=False, has_light_data=False,
              has_temporal_data=False):
        """
        Compute pitch

        Parameters
        ----------
        data_arr: np.ndarray [shape=(..., n)]
            Input audio array

        has_corr_data: bool
            If true, the `extra_data_dic` in the result will contain `corr_data` information.

            `corr_data` includes four arrays: `corr_arr`, `db_arr`, `height_arr`, and `len_arr`.

        has_cut_data: bool
            If true, the `extra_data_dic` in the result will contain `cut_data` information.
            `cut_data` contains only the first four sets of data from `corr_data`.

            `cut_data` includes four arrays: `corr_arr`, `db_arr`, `height_arr`, and `len_arr`.

        has_flag_data: bool
            If true, the `extra_data_dic` in the result will contain `flag_data` information.

            `flag_data` contains one array, `flag_arr`.

        has_light_data: bool
            If true, the `extra_data_dic` in the result will contain `light_data` information.

            `light_data` contains one array, `light_arr`.

        has_temporal_data: bool
            If true, the `extra_data_dic` in the result will contain `temporal_data` information.

            `temporal_data` includes three arrays: `avg_temp_arr`, `max_temp_arr`, `percent_temp_arr`.

        Returns
        -------
        fre_arr: np.ndarray [shape=(..., time)]
            frequency array. No data.

        db_arr: np.ndarray [shape=(..., time)]
            db array

        extra_data_dic: (optional) dict
            If all `has_xx_data` values are False, the `extra_data_dic` will not be included in the returned result.

            If any `has_xx_data` is True, the `extra_data_dic` will include the corresponding data.
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        has_format_data = False

        data_length = data_arr.shape[-1]
        corr_corr_arr, corr_db_arr, corr_height_arr, corr_len_arr = None, None, None, None
        cut_corr_arr, cut_db_arr, cut_height_arr, cut_len_arr = None, None, None, None
        flag_arr = None
        light_arr = None
        (f_flag_arr, f_fre_arr1, f_fre_arr2, f_fre_arr3,
         f_db_arr1, f_db_arr2, f_db_arr3) = None, None, None, None, None, None, None
        avg_temp_arr, max_temp_arr, percent_temp_arr = None, None, None
        if data_arr.ndim == 1:
            fre_arr, db_arr = self._pitch(data_arr)
            if has_corr_data:
                corr_corr_arr, corr_db_arr, corr_height_arr, corr_len_arr = self._get_corr_data(data_length)
            if has_cut_data:
                cut_corr_arr, cut_db_arr, cut_height_arr, cut_len_arr = self._get_cut_data(data_length)
            if has_flag_data:
                flag_arr = self._get_flag_data()
            if has_light_data:
                light_arr = self._get_light_data()
            if has_format_data:
                f_flag_arr, f_fre_arr1, f_fre_arr2, f_fre_arr3, f_db_arr1, f_db_arr2, f_db_arr3 = self._get_format_data(
                    data_length)
            if has_temporal_data:
                avg_temp_arr, max_temp_arr, percent_temp_arr = self._get_temporal_data()
        else:
            data_arr, o_channel_shape = format_channel(data_arr, 1)
            channel_num = data_arr.shape[0]

            fre_arr = []
            db_arr = []
            if has_corr_data:
                corr_corr_arr, corr_db_arr, corr_height_arr, corr_len_arr = [], [], [], []
            if has_cut_data:
                cut_corr_arr, cut_db_arr, cut_height_arr, cut_len_arr = [], [], [], []
            if has_flag_data:
                flag_arr = []
            if has_light_data:
                light_arr = []
            if has_format_data:
                (f_flag_arr, f_fre_arr1, f_fre_arr2, f_fre_arr3,
                 f_db_arr1, f_db_arr2, f_db_arr3) = [], [], [], [], [], [], []
            if has_temporal_data:
                avg_temp_arr, max_temp_arr, percent_temp_arr = [], [], []
            for i in range(channel_num):
                _fre_arr, _db_arr = self._pitch(data_arr[i])
                fre_arr.append(_fre_arr)
                db_arr.append(_db_arr)
                if has_corr_data:
                    _corr_corr_arr, _corr_db_arr, _corr_height_arr, _corr_len_arr = self._get_corr_data(data_length)
                    corr_corr_arr.append(_corr_corr_arr)
                    corr_db_arr.append(_corr_db_arr)
                    corr_height_arr.append(_corr_height_arr)
                    corr_len_arr.append(_corr_len_arr)
                if has_cut_data:
                    _cut_corr_arr, _cut_db_arr, _cut_height_arr, _cut_len_arr = self._get_cut_data(data_length)
                    cut_corr_arr.append(_cut_corr_arr)
                    cut_db_arr.append(_cut_db_arr)
                    cut_height_arr.append(_cut_height_arr)
                    cut_len_arr.append(_cut_len_arr)
                if has_flag_data:
                    _flag_arr = self._get_flag_data()
                    flag_arr.append(_flag_arr)
                if has_light_data:
                    _light_arr = self._get_light_data()
                    light_arr.append(_light_arr)
                if has_format_data:
                    (_f_flag_arr, _f_fre_arr1, _f_fre_arr2, _f_fre_arr3,
                     _f_db_arr1, _f_db_arr2, _f_db_arr3) = self._get_format_data(data_length)
                    f_flag_arr.append(_f_flag_arr)
                    f_fre_arr1.append(_f_fre_arr1)
                    f_fre_arr2.append(_f_fre_arr2)
                    f_fre_arr3.append(_f_fre_arr3)
                    f_db_arr1.append(_f_db_arr1)
                    f_db_arr2.append(_f_db_arr2)
                    f_db_arr3.append(_f_db_arr3)
                if has_temporal_data:
                    _avg_temp_arr, _max_temp_arr, _percent_temp_arr = self._get_temporal_data()
                    avg_temp_arr.append(_avg_temp_arr)
                    max_temp_arr.append(_max_temp_arr)
                    percent_temp_arr.append(_percent_temp_arr)

            fre_arr = np.stack(fre_arr, axis=0)
            db_arr = np.stack(db_arr, axis=0)
            fre_arr = revoke_channel(fre_arr, o_channel_shape, 1)
            db_arr = revoke_channel(db_arr, o_channel_shape, 1)
            if has_corr_data:
                corr_corr_arr = np.stack(corr_corr_arr, axis=0)
                corr_db_arr = np.stack(corr_db_arr, axis=0)
                corr_height_arr = np.stack(corr_height_arr, axis=0)
                corr_len_arr = np.stack(corr_len_arr, axis=0)

                corr_corr_arr = revoke_channel(corr_corr_arr, o_channel_shape, 1)
                corr_db_arr = revoke_channel(corr_db_arr, o_channel_shape, 1)
                corr_height_arr = revoke_channel(corr_height_arr, o_channel_shape, 1)
                corr_len_arr = revoke_channel(corr_len_arr, o_channel_shape, 1)

            if has_cut_data:
                cut_corr_arr = np.stack(cut_corr_arr, axis=0)
                cut_db_arr = np.stack(cut_db_arr, axis=0)
                cut_height_arr = np.stack(cut_height_arr, axis=0)
                cut_len_arr = np.stack(cut_len_arr, axis=0)

                cut_corr_arr = revoke_channel(cut_corr_arr, o_channel_shape, 1)
                cut_db_arr = revoke_channel(cut_db_arr, o_channel_shape, 1)
                cut_height_arr = revoke_channel(cut_height_arr, o_channel_shape, 1)
                cut_len_arr = revoke_channel(cut_len_arr, o_channel_shape, 1)

            if has_flag_data:
                flag_arr = np.stack(flag_arr, axis=0)

                flag_arr = revoke_channel(flag_arr, o_channel_shape, 1)

            if has_light_data:
                light_arr = np.stack(light_arr, axis=0)

                light_arr = revoke_channel(light_arr, o_channel_shape, 1)

            if has_format_data:
                f_flag_arr = np.stack(f_flag_arr, axis=0)
                f_fre_arr1 = np.stack(f_fre_arr1, axis=0)
                f_fre_arr2 = np.stack(f_fre_arr2, axis=0)
                f_fre_arr3 = np.stack(f_fre_arr3, axis=0)
                f_db_arr1 = np.stack(f_db_arr1, axis=0)
                f_db_arr2 = np.stack(f_db_arr2, axis=0)
                f_db_arr3 = np.stack(f_db_arr3, axis=0)

                f_flag_arr = revoke_channel(f_flag_arr, o_channel_shape, 1)
                f_fre_arr1 = revoke_channel(f_fre_arr1, o_channel_shape, 1)
                f_fre_arr2 = revoke_channel(f_fre_arr2, o_channel_shape, 1)
                f_fre_arr3 = revoke_channel(f_fre_arr3, o_channel_shape, 1)
                f_db_arr1 = revoke_channel(f_db_arr1, o_channel_shape, 1)
                f_db_arr2 = revoke_channel(f_db_arr2, o_channel_shape, 1)
                f_db_arr3 = revoke_channel(f_db_arr3, o_channel_shape, 1)

            if has_temporal_data:
                avg_temp_arr = np.stack(avg_temp_arr, axis=0)
                max_temp_arr = np.stack(max_temp_arr, axis=0)
                percent_temp_arr = np.stack(percent_temp_arr, axis=0)

                avg_temp_arr = revoke_channel(avg_temp_arr, o_channel_shape, 1)
                max_temp_arr = revoke_channel(max_temp_arr, o_channel_shape, 1)
                percent_temp_arr = revoke_channel(percent_temp_arr, o_channel_shape, 1)

        extra_data_dic = {}
        if has_corr_data:
            extra_data_dic['corr_data'] = (corr_corr_arr, corr_db_arr, corr_height_arr, corr_len_arr)
        if has_cut_data:
            extra_data_dic['cut_data'] = (cut_corr_arr, cut_db_arr, cut_height_arr, cut_len_arr)
        if has_flag_data:
            extra_data_dic['flag_data'] = (flag_arr,)
        if has_light_data:
            extra_data_dic['light_data'] = (light_arr,)
        if has_format_data:
            extra_data_dic['format_data'] = (f_flag_arr, f_fre_arr1, f_fre_arr2, f_fre_arr3,
                                             f_db_arr1, f_db_arr2, f_db_arr3)
        if has_temporal_data:
            extra_data_dic['temporal_data'] = (avg_temp_arr, max_temp_arr, percent_temp_arr)

        if extra_data_dic:
            return fre_arr, db_arr, extra_data_dic
        else:
            return fre_arr, db_arr

    # def debug(self, flag=True):
    #     fn = self._lib['pitchFFPObj_enableDebug']
    #     fn.argtypes = [POINTER(OpaquePitchFFP), c_int]
    #     fn(self._obj, c_int(int(flag)))

    def __del__(self):
        if self._is_created:
            fn = self._lib['pitchFFPObj_free']
            fn.argtypes = [POINTER(OpaquePitchFFP)]
            fn.restype = c_void_p
            fn(self._obj)
