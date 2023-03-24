import os
import warnings
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_char_p
import scipy
import numpy as np
import soundfile as sf

from audioflux.base import Base
from audioflux.utils.util import ascontiguous_T

__all__ = [
    'read',
    'write',
    'convert_mono',
    'resample',
    'chirp',
    'WaveReader',
    'WaveWriter'
]


def read(path=None, dir=None, is_mono=True, samplate=None, re_type='scipy'):
    """
    Load an audio file as a NumPy array of floats.

    Parameters
    ----------
    path: str or list
        Audio file path or path list.

        If `list`, the sampling rate and audio length of all audio files must be the same

    dir: str
        Audio file directory.

        If you set the directory, the `path` parameters cannot work.

        The sampling rate and audio length in the directory must be the same.

    is_mono: bool
        Whether to convert it into a single channel

    samplate: int or None
        Convert audio sampling rate. None defaults to use audio

    re_type: str
        Resample type

        - scipy: scipy.signal.resample
        - scipy_poly: scipy.signal.resample_poly

    Returns
    -------
    y: np.ndarray [shape=(..., n)]
        audio time series.

    sr: number > 0 [scalar]
        sampling rate of y
    """

    if dir is not None:
        path = []
        for filename in os.listdir(dir):
            path.append(os.path.join(dir, filename))

    if isinstance(path, str):
        data, sr = __read(path)
        if is_mono:
            data = convert_mono(data)
        elif data.ndim == 1:
            data = data.reshape(1, -1)
    else:
        sr = None
        data = []
        audio_shape = None
        for fp in path:
            try:
                _data, _sr = __read(fp)
            except Exception as e:
                warnings.warn(f'Load file error, skip: {fp}, {e}')
                continue

            # check sampling rate
            if sr is None:
                sr = _sr
            elif sr != _sr:
                raise ValueError('When loading multiple audio files, the sampling rate must be the same')

            # check shape
            if audio_shape is None:
                audio_shape = _data.shape
            elif audio_shape != _data.shape:
                raise ValueError('When loading multiple audio files, the audio shape must be the same')

            if is_mono:
                _data = convert_mono(_data)
            elif _data.ndim == 1:
                _data = _data.reshape(1, -1)

            data.append(_data)
        data = np.stack(data, axis=0)

    if samplate is not None and samplate != sr:
        data = resample(data, source_samplate=sr, target_samplate=samplate, re_type=re_type)
        sr = samplate
    return data, sr


def __read(fp):
    with sf.SoundFile(fp, 'r') as f:
        data = f.read()
        data = data.astype(dtype=np.float32)
        sr = f.samplerate
    data = ascontiguous_T(data)
    return data, sr


def write(path, data, samplate=32000, subtype='PCM_32', format='WAV'):
    """
    Write audio data from a NumPy array to the file.

    Parameters
    ----------
    path: str
        Path to save audio file.

    data: np.ndarray [shape=(frames,) or (channel, frames)]
        Audio file data

    samplate: int
        Audio file sample rate

    subtype: str
        The subtype of the sound file.

        See: `soundfile.SoundFile`
    format: str
        The major format of the sound file.

        See: `soundfile.SoundFile`
    """
    data = np.asarray(data, dtype=np.float32, order='C')

    if data.ndim > 2:
        raise ValueError(f'data must be less than equal to 2 dimensions')

    channel = 1 if data.ndim == 1 else data.shape[0]
    if channel != 1:
        data = data.T

    with sf.SoundFile(path, 'w', samplerate=samplate, channels=channel,
                      subtype=subtype, format=format) as f:
        f.write(data)


def convert_mono(x):
    """
    Convert audio data from multi-channel to single-channel

    Parameters
    ----------
    x: np.ndarray [shape=(channel, frames) or (samples, channel, frames)]
        Audio file data of multi-channel

    Returns
    -------
    out: np.ndarray [shape(=(frames,)]
        Audio file data of single-channel
    """

    if x.ndim > 1:
        x = np.mean(x, axis=-2)
    return x


def resample(x, source_samplate, target_samplate, re_type='scipy'):
    """
    Resample audio data from source_sr to target_sr

    Parameters
    ----------
    x: np.ndarray [shape=(..., frames)]
        Audio data

    source_samplate: int
        Audio's source sampling rate

    target_samplate: int
        Audio's target sampling rate

    re_type: str
        Resample type

        - scipy: scipy.signal.resample
        - scipy_poly: scipy.signal.resample_poly

    Returns
    -------
    out: np.ndarray [shape=(..., frames)]
        Resample audio data
    """
    x = np.asarray(x, dtype=np.float32, order='C')

    if target_samplate == source_samplate:
        return x

    if not 8000 <= target_samplate < source_samplate:
        raise ValueError(
            f'target_samplate[{target_samplate}] must be between 8000 to source_samplate[{source_samplate}]')

    if re_type == 'scipy':
        num = int(np.ceil(x.shape[-1] * (target_samplate * 1.0 / source_samplate)))
        y = scipy.signal.resample(x, num, axis=-1)
    elif re_type == 'scipy_poly':
        gcd = np.gcd(source_samplate, target_samplate)
        up = target_samplate // gcd
        down = source_samplate // gcd
        y = scipy.signal.resample_poly(x, up=up, down=down, axis=-1)
    else:
        raise ValueError(f're_type[{re_type}] not supported')
    return y


def chirp(fmin, fmax, duration, samplate=32000, phi=None, method='logarithmic'):
    """
    chirp signal

    The chirp sweeps from frequency ``fmin`` to ``fmax`` (in Hz).

    Parameters
    ----------
    fmin: float
        initial frequency

    fmax: float
        final frequency

    duration: float
        output duration in seconds.

    samplate: int
        desired sampling rate of the output signal

    phi: float, optional
        Phase offset, in radians.

        If None, use `1 / 2 * -np.pi`.

    method: str
        linear, quadratic, logarithmic, hyperbolic

        See: `scipy.signal.chirp`

    Returns
    -------
    out: np.ndarray [shape=(duration*samplate, )]
        chirp single array
    """
    if fmin <= 0 or fmax <= 0:
        raise ValueError(f'fmax and fmin must be strictly positive')

    t = np.arange(duration, step=1. / samplate)
    if phi is None:
        phi = 1 / 2 * -np.pi
    y = scipy.signal.chirp(t, fmin, duration, fmax, method=method, phi=phi / np.pi * 180)

    return y


class OpaqueWaveRead(Structure):
    _fields_ = []


class OpaqueWaveWrite(Structure):
    _fields_ = []


class WaveReader(Base):
    """
    Wave Reader

    Parameters
    ----------
    file_path: str
        file path
    """

    def __init__(self, file_path):
        super(WaveReader, self).__init__(pointer(OpaqueWaveRead()))
        self.file_path = file_path

        fn = self._lib['waveReadObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueWaveRead)), c_char_p]
        fn(self._obj, self.file_path.encode())
        self._is_created = True

    def get_infor(self) -> dict:
        """
        Get wave info

        Returns
        -------
        out: dict
            samplate: int
            bit: int
            channel_num: int
        """
        fn = self._lib['waveReadObj_getInfor']
        fn.argtypes = [
            POINTER(OpaqueWaveRead),
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_int),
        ]
        samplate_c = c_int(0)
        bit_c = c_int(0)
        channel_num_c = c_int(0)
        fn(self._obj,
           pointer(samplate_c),
           pointer(bit_c),
           pointer(channel_num_c))
        return {'samplate': samplate_c.value,
                'bit': bit_c.value,
                'channel_num': channel_num_c.value}

    def read(self, n):
        """
        Read at most n characters from stream.

        Parameters
        ----------
        n: int

        Returns
        -------
        out: np.ndarray [shape=(data,) or (channel, data)]
        """
        channel_num = self.get_infor()['channel_num']

        fn = self._lib['waveReadObj_read']
        fn.argtypes = [
            POINTER(OpaqueWaveRead),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int,
        ]
        _data_len = n * channel_num
        audio_arr = np.zeros(_data_len, dtype=np.float32)
        fn(self._obj, audio_arr, _data_len)
        if channel_num != 1:
            audio_arr = audio_arr.reshape((channel_num, n), order='F')
            audio_arr = np.ascontiguousarray(audio_arr)
        return audio_arr

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['waveReadObj_free']
            free_fn.argtypes = [POINTER(OpaqueWaveRead)]
            free_fn.restype = c_void_p
            free_fn(self._obj)


class WaveWriter(Base):
    """
    Wave Writer

    Parameters
    ----------
    file_path: str
        Wave writer path

    samplate: int
        Audio file sample rate

    bit: int
        audio bit depth

    channel_num: int
        Number of audio channels. 1 is mono, 2 is stereo
    """

    def __init__(self, file_path, samplate=32000, bit=16, channel_num=1):
        super(WaveWriter, self).__init__(pointer(OpaqueWaveWrite()))
        self.file_path = file_path
        self.samplate = samplate
        self.bit = bit
        self.channel_num = channel_num

        fn = self._lib['waveWriteObj_new']
        fn.argtypes = [POINTER(POINTER(OpaqueWaveWrite)),
                       c_char_p,
                       POINTER(c_int),
                       POINTER(c_int),
                       POINTER(c_int)]
        fn(self._obj,
           self.file_path.encode(),
           pointer(c_int(samplate)),
           pointer(c_int(bit)),
           pointer(c_int(channel_num)))
        self._is_created = True

    def write(self, data_arr):
        """
        Write data array to a wave file.

        Parameters
        ----------
        data_arr: np.ndarray [shape=(data,) or (channel, data)]

        Returns
        -------

        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        channel = 0
        data_len = 0
        if data_arr.ndim == 1:
            channel = 1
            data_len = data_arr.shape[-1]
        elif data_arr.ndim == 2:
            channel = data_arr.shape[-2]
            data_len = data_arr.shape[-1]

        if channel == 0 or channel != self.channel_num:
            raise ValueError(f"Invalid shape: {data_arr.shape}")

        data_arr = data_arr.reshape(channel * data_len, order='F')
        data_arr = np.ascontiguousarray(data_arr, dtype=np.float32)

        fn = self._lib['waveWriteObj_write']
        data_arr = data_arr.astype(np.float32)
        fn.argtypes = [
            POINTER(OpaqueWaveWrite),
            np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
            c_int,
        ]
        fn(self._obj, data_arr, data_arr.size)

    def __del__(self):
        if self._is_created:
            free_fn = self._lib['waveWriteObj_free']
            free_fn.argtypes = [POINTER(OpaqueWaveWrite)]
            free_fn.restype = c_void_p
            free_fn(self._obj)
