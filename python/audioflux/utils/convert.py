from ctypes import c_int, c_float, POINTER, pointer
from functools import reduce
import numpy as np

from audioflux.fftlib import get_fft_lib
from audioflux.utils.util import ascontiguous_swapaxex

__all__ = [
    'power_to_db',
    'power_to_abs_db',
    'mag_to_abs_db',
    'log_compress',
    'log10_compress',
    'delta',
    'get_phase',
    'note_to_midi',
    'midi_to_hz',
    'note_to_hz',
    'midi_to_note',
    'hz_to_midi',
    'hz_to_note',
]


def power_to_db(X, min_db=-80):
    """
    Convert power spectrogram to relative decibel(dB) spectrogram.

    Parameters
    ----------
    X: np.ndarray [shape=(..., fre, time)]
        Input array

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        relative dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    if X.ndim < 2:
        raise ValueError('The dimension should be greater than equal to 2')
    fn = get_fft_lib()['util_powerToDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    if X.ndim == 2:
        shape = X.shape
        arr = X.flatten()
        length = arr.size
        db_arr = np.zeros(arr.shape, dtype=np.float32)
        fn(arr, c_int(length), c_float(min_db), db_arr)
        db_arr = db_arr.reshape(shape)
    else:
        *channel_shape, n_fre, n_time = X.shape
        channel_num = reduce(lambda x, y: x * y, channel_shape)
        arr = X.reshape(channel_num, n_fre * n_time)
        db_arr = np.zeros_like(arr, dtype=arr.dtype)
        length = arr.shape[-1]
        for i in range(channel_num):
            fn(arr[i], c_int(length), c_float(min_db), db_arr[i])
        db_arr = db_arr.reshape((*channel_shape, n_fre, n_time))
    return db_arr


def power_to_abs_db(X, fft_length=4096, is_norm=False, min_db=-80):
    """
    Convert power spectrogram to absolute decibel(dB) spectrogram

    Parameters
    ----------
    X: np.ndarray [shape=(..., fre, time)]
        Input array

    fft_length: int
        fft length

    is_norm: bool
        Whether to use normalization

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        absolute dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_powerToAbsDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    length = arr.size
    db_arr = np.zeros(arr.shape, dtype=np.float32)
    fn(arr, c_int(length), c_int(fft_length), c_int(int(is_norm)), c_float(min_db), db_arr)
    return db_arr.reshape(shape)


def mag_to_abs_db(X, fft_length=4096, is_norm=False, min_db=-80):
    """
    Convert magnitude spectrogram to absolute decibel(dB) spectrogram

    Parameters
    ----------
    X: np.ndarray [shape=(..., fre, time)]
        Input array

    fft_length: int
        fft length

    is_norm: bool
        Whether to use normalization

    min_db: float
        Minimum dB. Values smaller than `min_db` will be set to `min_db`

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        absolute dB array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    fn = get_fft_lib()['util_magToAbsDB']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        c_int,
        c_float,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    length = arr.size
    db_arr = np.zeros(arr.shape, dtype=np.float32)
    fn(arr, c_int(length), c_int(fft_length), c_int(int(is_norm)), c_float(min_db), db_arr)
    return db_arr.reshape(shape)


def log_compress(X, gamma=1.0):
    """
    log compression

    Parameters
    ----------
    X: np.ndarray [shape=(..., fre,) or (..., fre, time)]
        Input array

    gamma: float
        1/10/20/...

    Returns
    -------
    out: np.ndarray [shape=(..., fre,) or (..., fre, time)]
        log compressed array
    """
    X = np.asarray(X, dtype=np.float32, order='C')

    fn = get_fft_lib()['util_logCompress']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        POINTER(c_float),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    ret_arr = np.zeros_like(arr, dtype=np.float32)
    fn(arr, pointer(c_float(gamma)), c_int(arr.shape[-1]), ret_arr)
    ret_arr = ret_arr.reshape(shape)
    return ret_arr


def log10_compress(X, gamma=1.0):
    """
    log10 compression

    Parameters
    ----------
    X: np.ndarray [shape=(..., fre,) or (..., fre, time)]
        Input array

    gamma: float
        1/10/20/...

    Returns
    -------
    out: np.ndarray [shape=(..., fre,) or (..., fre, time)]
        log10 compressed array
    """
    X = np.asarray(X, dtype=np.float32, order='C')

    fn = get_fft_lib()['util_log10Compress']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        POINTER(c_float),
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    arr = X.flatten()
    ret_arr = np.zeros_like(arr, dtype=np.float32)
    fn(arr, pointer(c_float(gamma)), c_int(arr.shape[-1]), ret_arr)
    ret_arr = ret_arr.reshape(shape)
    return ret_arr


def delta(X, order=9):
    """
    Compute delta features

    Parameters
    ----------
    X: np.ndarray [shape=(..., n_fre, n_time)]
        Input array

    order: int
        must odd

    Returns
    -------
    out: np.ndarray [shape=(..., n_fre, n_time)]
        delta array
    """
    X = np.asarray(X, dtype=np.float32, order='C')
    if X.ndim < 2:
        raise ValueError('The dimension should be greater than equal to 2')

    X = ascontiguous_swapaxex(X, -1, -2)
    fn = get_fft_lib()['util_delta']
    fn.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
        c_int,
        c_int,
        np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
    ]

    shape = X.shape
    *num, n_fre = X.shape
    num = reduce(lambda x, y: x * y, num)
    arr = X.reshape((num, n_fre))
    ret_arr = np.zeros_like(arr, dtype=np.float32)
    for i in range(num):
        fn(arr[i], c_int(n_fre), c_int(order), ret_arr[i])
    ret_arr = ret_arr.reshape(shape)
    ret_arr = ascontiguous_swapaxex(ret_arr, -1, -2)
    return ret_arr


def get_phase(D):
    """
    Extract phase data from a complex-valued spectrogram D.

    Parameters
    ----------
    D: np.ndarray [shape=(..., fre, time)]
        a complex-valued spectrogram

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        phase data
    """
    if not np.iscomplexobj(D):
        raise ValueError(f'D must be type of complex')

    m_real_arr = D.real.copy()
    m_imag_arr = D.imag.copy()
    m_real_arr[m_real_arr < 1e-16] = 1e-16
    m_phase_arr = np.arctan2(m_imag_arr, m_real_arr)
    return m_phase_arr


def note_to_midi(note):
    """
    Convert note to MIDI number.

    Parameters
    ----------
    note: str
        note names.

    Returns
    -------
    out: int
        MIDI number

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.note_to_midi('C1')
    24
    >>> af.utils.note_to_midi('D#2')
    39
    >>> af.utils.note_to_midi('#F4')
    66
    >>> af.utils.note_to_midi('Gb4')
    66
    """
    if not isinstance(note, str):
        raise ValueError(f'note must be type of str')

    pitch_dic = {'C': 0,
                 'D': 2,
                 'E': 4,
                 'F': 5,
                 'G': 7,
                 'A': 9,
                 'B': 11}
    acc_dic = {
        '#': 1,
        'b': -1
    }

    pitch = None
    acc = 0
    octave = 0
    for s in note:
        if s in pitch_dic:
            if pitch is not None:
                raise ValueError(f'note={note} formatting failure')
            pitch = pitch_dic[s]
        elif s in acc_dic:
            acc += acc_dic[s]
        elif s.isdigit():
            octave = octave * 10 + int(s)
        else:
            raise ValueError(f'note={note} formatting failure')
    if pitch is None:
        raise ValueError(f'note={note} formatting failure')

    midi = pitch + acc + 12 * octave + 12
    return midi


def midi_to_hz(midi):
    """
    Convert MIDI number to the frequency (Hz).

    Parameters
    ----------
    midi: int
        MIDI number.

    Returns
    -------
    out: float
        frequency (Hz)

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.midi_to_hz(24)
    32.70319566257483
    >>> af.utils.midi_to_hz(69)
    440.0
    >>> af.utils.midi_to_hz(81)
    880.0
    """
    return 440.0 * (2.0 ** ((np.asanyarray(midi) - 69.0) / 12.0))


def note_to_hz(note):
    """
    Convert note to the frequency (Hz).

    Parameters
    ----------
    note: str
        note name

    Returns
    -------
    out: float
        frequency (Hz)

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.note_to_hz('C1')
    32.70319566257483
    >>> af.utils.note_to_hz('G#3')
    207.65234878997256
    >>> af.utils.note_to_hz('A4')
    440.0
    >>> af.utils.note_to_hz('Gb5')
    739.9888454232688

    """
    return midi_to_hz(note_to_midi(note))


def midi_to_note(midi, is_octave=True):
    """
    Convert MIDI number to note.

    Parameters
    ----------
    midi: int or float
        MIDI number

    is_octave: bool
        If `Ture`, then show the octave number

    Returns
    -------
    out: str
        note name

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.midi_to_note(0)
    'C-1'
    >>> af.utils.midi_to_note(24)
    'C1'
    >>> af.utils.midi_to_note(69.4)
    'A4'
    >>> af.utils.midi_to_note(69.7)
    'A#4'
    >>> af.utils.midi_to_note(90, octave=False)
    'F#'
    """
    midi = round(midi)
    pitch_list = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    octave = midi // 12 - 1
    pitch = midi % 12

    if is_octave:
        note = f'{pitch_list[pitch]}{octave}'
    else:
        note = f'{pitch_list[pitch]}'

    return note


def hz_to_midi(frequencies):
    """
    Convert the frequency (Hz) to MIDI number.

    Parameters
    ----------
    frequencies: float
        the frequency

    Returns
    -------
    out:
        MIDI number

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.hz_to_midi(32.703)
    23.99989642029834
    >>> af.utils.hz_to_midi(110)
    45.0
    >>> af.utils.hz_to_midi(220)
    57.0
    >>> af.utils.hz_to_midi(440)
    69.0
    """
    return 12 * (np.log2(np.asanyarray(frequencies)) - np.log2(440.0)) + 69


def hz_to_note(frequencies):
    """
    Convert the frequency (Hz) to note.

    Parameters
    ----------
    frequencies: float
        the frequency

    Returns
    -------
    out:
        note

    Examples
    --------
    >>> import audioflux as af
    >>> af.utils.hz_to_note(32.703)
    'C1'
    >>> af.utils.hz_to_note(440)
    'A4'
    """
    midi = hz_to_midi(frequencies)
    note = midi_to_note(midi)
    return note
