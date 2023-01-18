import warnings
import numpy as np

__all__ = [
    'note_to_midi',
    'midi_to_hz',
    'note_to_hz',
    'midi_to_note',
    'hz_to_midi',
    'hz_to_note',
    'ascontiguous_T',
    'check_audio',
    'check_audio_length'
]


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
    >>> audioFlux.utils.note_to_midi('C1')
    24
    >>> audioFlux.utils.note_to_midi('D#2')
    39
    >>> audioFlux.utils.note_to_midi('#F4')
    66
    >>> audioFlux.utils.note_to_midi('Gb4')
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
    >>> audioflux.utils.midi_to_hz(24)
    32.70319566257483
    >>> audioflux.utils.midi_to_hz(69)
    440.0
    >>> audioflux.utils.midi_to_hz(81)
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
    >>> audioflux.utils.note_to_hz('C1')
    32.70319566257483
    >>> audioflux.utils.note_to_hz('G#3')
    207.65234878997256
    >>> audioflux.utils.note_to_hz('A4')
    440.0
    >>> audioflux.utils.note_to_hz('Gb5')
    739.9888454232688

    """
    return midi_to_hz(note_to_midi(note))


def midi_to_note(midi, is_octave=True):
    """
    Convert MIDI number to note.

    Parameters
    ----------
    midi: int
        MIDI number

    is_octave: bool
        If `Ture`, then show the octave number

    Returns
    -------
    out: str
        note name

    Examples
    --------
    >>> auidoflux.utils.midi_to_note(0)
    'C-1'
    >>> auidoflux.utils.midi_to_note(24)
    'C1'
    >>> auidoflux.utils.midi_to_note(69.4)
    'A4'
    >>> auidoflux.utils.midi_to_note(69.7)
    'A#4'
    >>> auidoflux.utils.midi_to_note(90, octave=False)
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
    >>> audioflux.utils.hz_to_midi(32.703)
    23.99989642029834
    >>> audioflux.utils.hz_to_midi(110)
    45.0
    >>> audioflux.utils.hz_to_midi(220)
    57.0
    >>> audioflux.utils.hz_to_midi(440)
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
    >>> audioflux.utils.hz_to_note(32.703)
    'C1'
    >>> audioflux.utils.hz_to_note(440)
    'A4'
    """
    midi = hz_to_midi(frequencies)
    note = midi_to_note(midi)
    return note


def ascontiguous_T(X, dtype=None, *args, **kwargs):
    """
    The transposed array, and then return a contiguous array in memory.

    Parameters
    ----------
    X: ndarray
        Input array.
    dtype: str or dtype object, optional
        Data-type of returned array.

    Returns
    -------
    out: ndarray
        return a transposed and contiguous array
    """
    return np.ascontiguousarray(X.T, dtype=dtype, *args, **kwargs)


def check_audio(X, is_mono=True):
    """
    check audio data

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio array

    is_mono: bool
        is mono

    Returns
    -------
    out: bool
    """
    if not isinstance(X, np.ndarray):
        raise ValueError(f'Audio data must be of type np.ndarray')

    if not X.flags.c_contiguous:
        raise ValueError("Audio data must be C-contiguous")

    if X.ndim == 0:
        raise ValueError(f'Audio data must have at least one dimension')

    if is_mono and X.ndim != 1:
        raise ValueError(f'Audio data must be a monophonic, there is only 1 dimension'
                         f', given X.ndim={X.ndim}, X.shape={X.shape}')

    return True


def check_audio_length(X, radix2_exp):
    data_len = X.shape[0]
    fft_length = 1 << radix2_exp
    if data_len < fft_length:
        pad_len = fft_length - data_len
        warnings.warn(f'The audio length={data_len} is not enough for fft_length={fft_length}(2**radix2_exp), '
                      f'and {pad_len} zeros are automatically filled after the audio')
        X = np.pad(X, (0, pad_len))
    elif data_len > fft_length:
        warnings.warn(f'fft_length={fft_length}(2**radix2_exp) is too small for data_arr length={data_len}, '
                      f'only the first fft_length={fft_length} data are valid')
        X = X[:fft_length].copy()
    return X
