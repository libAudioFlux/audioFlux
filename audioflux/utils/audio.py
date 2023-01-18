import numpy as np
from ctypes import Structure, POINTER, pointer, c_int, c_void_p, c_char_p
from audioflux.base import Base
from soundfile import SoundFile

__all__ = [
    'read',
    'write',
    'chirp',
    'WaveReader',
    'WaveWriter'
]


def read(path):
    """
    Load an audio file as a NumPy array of floats.

    Parameters
    ----------
    path: str
        audio file path

    Returns
    -------
    y: np.ndarray [shape=(n,)]
        audio time series.

    sr: number > 0 [scalar]
        sampling rate of y
    """

    with SoundFile(path, 'r') as f:
        data = f.read(dtype=np.float32)
    return data, f.samplerate


def write(path, data, samplate=32000, subtype='PCM_32', format='WAV'):
    """
    Write audio data from a NumPy array to the file.

    Parameters
    ----------
    path: str
        Path to save audio file.

    data: np.ndarray [shape=(channel, frames)]
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

    channel = 1 if data.ndim == 1 else data.shape[0]
    if channel != 1:
        data = data.T

    with SoundFile(path, 'w', samplerate=samplate, channels=channel,
                   subtype=subtype, format=format) as f:
        f.write(data)


def chirp(fmin, fmax, duration, samplate=32000, phi=0, linear=False):
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

    linear: boolean
        - If is True,  use ``f(t) = f0 + (f1 - f0) * t / t1``
        - If is False, use ``f(t) = f0 * (f1/f0)**(t/t1)``

    Returns
    -------
    out: np.ndarray [shape=(duration*samplate, )]
        chirp single array
    """
    if fmin <= 0 or fmax <= 0:
        raise ValueError(f'fmax and fmin must be strictly positive')
    t = np.arange(duration, step=1. / samplate)

    if linear:
        k = (fmax - fmin) / duration
        phase = 2 * np.pi * (fmin * t + 0.5 * k * t * t)
    else:
        if fmin == fmax:
            phase = 2 * np.pi * fmin * t
        else:
            beta = duration / np.log(fmax / fmin)
            phase = 2 * np.pi * beta * fmin * (np.power(fmax / fmin, t / duration) - 1.0)

    phi *= np.pi / 180

    return np.cos(phase + phi)


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
        fn(self._obj,
               self.file_path.encode())
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
            data_len = data_arr.shape[0]
        elif data_arr.ndim == 2:
            channel = data_arr.shape[0]
            data_len = data_arr.shape[1]

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
