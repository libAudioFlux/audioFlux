import os

__all__ = ['sample_path']

SAMPLE_DATA_PATH = os.path.join(os.path.abspath(__file__).rsplit(os.path.sep, 1)[0],
                                'sample_data')


def sample_path(name):
    """
    Get sample path

    Parameters
    ----------
    name: str
        sample name

        - '220' - guitar 220hz's audio file
        - '880' - guitar 880hz's audio file
        - 'guitar_chord1' - Audio files of guitar chords
        - 'guitar_chord2' - Audio files of guitar chords
        - 'chord_metronome1' - Audio file of guitar chords with metronome
        - 'chord_metronome2' - Audio file of guitar chords with metronome
        - 'voice' - Vocal Scale audio file

    Returns
    -------
    out: str
        audio file abspath
    """
    return os.path.join(SAMPLE_DATA_PATH, f'{name}.wav')
