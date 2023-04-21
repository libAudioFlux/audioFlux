import numpy as np
from audioflux import BFT, XXCC, CQT
from audioflux.spectrogram import Spectrogram, MelSpectrogram, BarkSpectrogram, ErbSpectrogram
from audioflux.type import WindowType, SpectralFilterBankScaleType, SpectralFilterBankStyleType, \
    SpectralFilterBankNormalType, SpectralDataType, SpectralFilterBankType, CepstralRectifyType
from audioflux.utils import note_to_hz

__all__ = ['linear_spectrogram',
           'mel_spectrogram',
           'bark_spectrogram',
           'erb_spectrogram',
           'cqt', 'vqt',
           'mfcc', 'bfcc', 'gtcc', 'cqcc',
           'chroma_linear', 'chroma_octave', 'chroma_cqt']


def linear_spectrogram(X, num=None, radix2_exp=12, samplate=32000, slide_length=None,
                       low_fre=0., window_type=WindowType.HANN,
                       style_type=SpectralFilterBankStyleType.SLANEY,
                       data_type=SpectralDataType.POWER,
                       is_reassign=False):
    """
    Short-time Fourier transform (Linear/STFT)

    It is a Fourier-related transform used to determine
    the sinusoidal frequency and phase content of
    local sections of a signal as it changes over time.

    .. Note:: We recommend using the `BFT` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int or None
        Number of frequency bins to generate, starting at `low_fre`.

        If `num` is None, then ``num = fft_length / 2 + 1``.
        When `radix2_exp=12, fft_length=4096, samplate=32000, low_fre=0`,
        then `num=2049`, and frequency range is 0-16000.

        The size of each band is ``samplate / fft_length``.

    radix2_exp: int
        ``fft_length=2**radix2_exp``.

    samplate: int
        Sampling rate of the incoming audio.

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    low_fre: float
        Lowest frequency.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    is_reassign: bool
        Whether to use reassign.

    Returns
    -------
    spectrogram: np.ndarray [shape=(..., fre, time)]
        The matrix of Linear(STFT)

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    mel_spectrogram
    bark_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = 0
    >>> spec_arr, fre_band_arr = af.linear_spectrogram(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Linear Spectrogram')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")
    """
    if num is None:
        num = (1 << radix2_exp) // 2 + 1
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.LINEAR, style_type=style_type,
                  data_type=data_type, is_reassign=is_reassign)
    spec_arr = bft_obj.bft(X, result_type=1)
    fre_band_arr = bft_obj.get_fre_band_arr()
    return spec_arr, fre_band_arr


def mel_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                    low_fre=0., high_fre=None,
                    window_type=WindowType.HANN, slide_length=None,
                    style_type=SpectralFilterBankStyleType.SLANEY,
                    normal_type=SpectralFilterBankNormalType.NONE,
                    data_type=SpectralDataType.POWER):
    """
    Mel-scale spectrogram.

    .. Note:: We recommend using the `BFT` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int
        Number of mel frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``.

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000 (samplate / 2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    Returns
    -------
    spectrogram: np.ndarray [shape=(..., fre, time)]
        The matrix of MEL

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    linear_spectrogram
    bark_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = 0
    >>> spec_arr, fre_band_arr = af.mel_spectrogram(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Mel Spectrogram')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")

    """
    spec_obj = MelSpectrogram(num=num, radix2_exp=radix2_exp, samplate=samplate,
                              low_fre=low_fre, high_fre=high_fre,
                              window_type=window_type, slide_length=slide_length,
                              style_type=style_type,
                              normal_type=normal_type, data_type=data_type)
    spec_arr = spec_obj.spectrogram(X)
    fre_band_arr = spec_obj.get_fre_band_arr()
    return spec_arr, fre_band_arr


def bark_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                     low_fre=None, high_fre=None,
                     window_type=WindowType.HANN, slide_length=None,
                     style_type=SpectralFilterBankStyleType.SLANEY,
                     normal_type=SpectralFilterBankNormalType.NONE,
                     data_type=SpectralDataType.POWER):
    """
    Bark-scale spectrogram.

    .. Note:: We recommend using the `BFT` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int
        Number of bark frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``.

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000 (samplate / 2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    Returns
    -------
    spectrogram: np.ndarray [shape=(..., fre, time)]
        The matrix of BARK

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    linear_spectrogram
    mel_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = 0
    >>> spec_arr, fre_band_arr = af.bark_spectrogram(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Bark Spectrogram')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")
    """
    spec_obj = BarkSpectrogram(num=num, radix2_exp=radix2_exp, samplate=samplate,
                               low_fre=low_fre, high_fre=high_fre,
                               window_type=window_type, slide_length=slide_length,
                               style_type=style_type,
                               normal_type=normal_type, data_type=data_type)
    spec_arr = spec_obj.spectrogram(X)
    fre_band_arr = spec_obj.get_fre_band_arr()
    return spec_arr, fre_band_arr


def erb_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                    low_fre=None, high_fre=None,
                    window_type=WindowType.HANN, slide_length=None,
                    style_type=SpectralFilterBankStyleType.SLANEY,
                    normal_type=SpectralFilterBankNormalType.NONE,
                    data_type=SpectralDataType.POWER):
    """
    Erb-scale spectrogram.

    .. Note:: We recommend using the `BFT` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int
        Number of erb frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``.

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000 (samplate / 2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

        See: `type.SpectralDataType`

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        The matrix of ERB

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    linear_spectrogram
    mel_spectrogram
    bark_spectrogram

    BFT
    NSGT
    CWT
    PWT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = 0
    >>> spec_arr, fre_band_arr = af.erb_spectrogram(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='Erb Spectrogram')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")
    """
    spec_obj = ErbSpectrogram(num=num, radix2_exp=radix2_exp, samplate=samplate,
                              low_fre=low_fre, high_fre=high_fre,
                              window_type=window_type, slide_length=slide_length,
                              style_type=style_type,
                              normal_type=normal_type, data_type=data_type)
    spec_arr = spec_obj.spectrogram(X)
    fre_band_arr = spec_obj.get_fre_band_arr()
    return spec_arr, fre_band_arr


def mfcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         mel_num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Mel-frequency cepstral coefficients (MFCCs)

    .. Note:: We recommend using the `BFT` and `XXCC` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    cc_num: int
        number of MFCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    mel_num: int
        Number of mel frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    low_fre: float or None
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Returns
    -------
    out: np.ndarray [shape=(..., cc_num, time)]
        The matrix of MFCCs

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of Mel frequency bands

    See Also
    --------
    bfcc
    gtcc
    cqcc

    BFT
    XXCC

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract mfcc data

    >>> cc_arr, _ = af.mfcc(audio_arr, samplate=sr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, cc_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cc_arr, axes=ax,
    >>>                 x_coords=x_coords, x_axis='time',
    >>>                 title='MFCC')
    >>> fig.colorbar(img, ax=ax)
    """
    bft_obj = BFT(num=mel_num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.MEL,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[-1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    fre_band_arr = bft_obj.get_fre_band_arr()
    return xxcc_arr, fre_band_arr


def bfcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         bark_num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Bark-frequency cepstral coefficients (BFCCs)

    .. Note:: We recommend using the `BFT` and `XXCC` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    cc_num: int
        number of BFCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    bark_num: int
        Number of bark frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    low_fre: float or None
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Returns
    -------
    out: np.ndarray [shape=(..., cc_num, time)]
        The matrix of BFCCs

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of Bark frequency bands

    See Also
    --------
    mfcc
    gtcc
    cqcc

    BFT
    XXCC

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract bfcc data

    >>> cc_arr, _ = af.bfcc(audio_arr, samplate=sr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, cc_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cc_arr, axes=ax,
    >>>                 x_coords=x_coords, x_axis='time',
    >>>                 title='BFCC')
    >>> fig.colorbar(img, ax=ax)
    """
    bft_obj = BFT(num=bark_num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.BARK,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[-1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    fre_band_arr = bft_obj.get_fre_band_arr()
    return xxcc_arr, fre_band_arr


def gtcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         erb_num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Gammatone cepstral coefficients (GTCCs)

    .. Note:: We recommend using the `BFT` and `XXCC` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    cc_num: int
        number of GTCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    erb_num: int
        Number of erb frequency bins to generate, starting at `low_fre`.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    low_fre: float or None
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    Returns
    -------
    out: np.ndarray [shape=(..., cc_num, time)]
        The matrix of GTCCs

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of Erb frequency bands

    See Also
    --------
    mfcc
    bfcc
    cqcc

    BFT
    XXCC

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract gtcc data

    >>> cc_arr, _ = af.gtcc(audio_arr, samplate=sr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, cc_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cc_arr, axes=ax,
    >>>                 x_coords=x_coords, x_axis='time',
    >>>                 title='GTCC')
    >>> fig.colorbar(img, ax=ax)
    """
    bft_obj = BFT(num=erb_num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.ERB,
                  style_type=SpectralFilterBankStyleType.GAMMATONE,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[-1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    fre_band_arr = bft_obj.get_fre_band_arr()
    return xxcc_arr, fre_band_arr


def cqcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         cqt_num=84, samplate=32000, low_fre=note_to_hz('C1'), slide_length=None,
         bin_per_octave=12, window_type=WindowType.HANN,
         normal_type=SpectralFilterBankNormalType.AREA,
         is_scale=True):
    """
    Constant-Q cepstral coefficients (CQCCs)

    .. Note:: We recommend using the :ref:`CQT <transforms/cqt:CQT>` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    cc_num: int
        number of GTCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    cqt_num: int
        Number of cqt frequency bins to generate, starting at `low_fre`.

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float or None
        Lowest frequency.

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    bin_per_octave: int
        Number of bins per octave.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(..., cc_num, time)]
        The matrix of CQCCs

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of cqt frequency bands

    See Also
    --------
    mfcc
    bfcc
    gtcc

    CQT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract cqcc data

    >>> cc_arr, _ = af.cqcc(audio_arr, samplate=sr)

    Show plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, cc_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(cc_arr, axes=ax,
    >>>                 x_coords=x_coords, x_axis='time',
    >>>                 title='CQCC')
    >>> fig.colorbar(img, ax=ax)
    """
    cqt_obj = CQT(num=cqt_num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale)
    spec_arr = cqt_obj.cqt(X)
    power_arr = np.abs(spec_arr) ** 2
    cqcc_arr = cqt_obj.cqcc(power_arr, cc_num=cc_num, rectify_type=rectify_type)
    fre_band_arr = cqt_obj.get_fre_band_arr()
    return cqcc_arr, fre_band_arr


def cqt(X, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
        factor=1., thresh=0.01,
        window_type=WindowType.HANN, slide_length=None,
        normal_type=SpectralFilterBankNormalType.AREA,
        is_scale=True):
    """
    Constant-Q transform (CQT)

    .. Note:: We recommend using the :ref:`CQT <transforms/cqt:CQT>` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int
        Number of frequency bins to generate, starting at `low_fre`.

        Usually: ``num = octave * bin_per_octave``, `default: 84 (7 * 12)`

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    bin_per_octave: int
        Number of bins per octave.

    factor: float
        Factor value

    thresh: float
        Thresh value

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        The matrix of CQT

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    cqcc
    vqt

    CQT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = af.utils.note_to_hz('C1')
    >>> spec_arr, fre_band_arr = af.cqt(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr ** 2)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='CQT')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=0., thresh=thresh)
    cqt_arr = cqt_obj.cqt(X)
    cqt_arr = np.abs(cqt_arr)
    fre_band_arr = cqt_obj.get_fre_band_arr()
    return cqt_arr, fre_band_arr


def vqt(X, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
        factor=1., beta=0.5, thresh=0.01,
        window_type=WindowType.HANN, slide_length=None,
        normal_type=SpectralFilterBankNormalType.AREA,
        is_scale=True):
    """
    Variable-Q transform (VQT)

    .. Note:: We recommend using the :ref:`CQT <transforms/cqt:CQT>` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    num: int
        Number of frequency bins to generate, starting at `low_fre`.

        Usually: ``num = octave * bin_per_octave``, `default: 84 (7 * 12)`

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    bin_per_octave: int
        Number of bins per octave.

    factor: float
        Factor value

    beta: float
        Beta value

    thresh: float
        Thresh value

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(..., fre, time)]
        The matrix of VQT

    fre_band_arr: np:ndarray [shape=(fre,)]
        The array of frequency bands

    See Also
    --------
    cqcc
    vqt

    CQT

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract spectrogram of dB

    >>> low_fre = af.utils.note_to_hz('C1')
    >>> spec_arr, fre_band_arr = af.vqt(audio_arr, samplate=sr, low_fre=low_fre)
    >>> spec_dB_arr = af.utils.power_to_db(spec_arr ** 2)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x/y-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, spec_arr.shape[-1] + 1)
    >>> y_coords = np.insert(fre_band_arr, 0, low_fre)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(spec_dB_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 y_coords=y_coords,
    >>>                 x_axis='time', y_axis='log',
    >>>                 title='VQT')
    >>> fig.colorbar(img, ax=ax, format="%+2.0f dB")
    """
    vqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=beta, thresh=thresh)
    vqt_arr = vqt_obj.cqt(X)
    vqt_arr = np.abs(vqt_arr)
    fre_band_arr = vqt_obj.get_fre_band_arr()
    return vqt_arr, fre_band_arr


def chroma_linear(X, chroma_num=12, radix2_exp=12, samplate=32000, low_fre=0., high_fre=16000.,
                  window_type=WindowType.HANN, slide_length=1024, data_type=SpectralDataType.POWER,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.NONE):
    """
    Linear(STFT) chromagram

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    chroma_num: int
        Number of chroma bins to generate.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency.

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    data_type: SpectralDataType
        Spectrogram data type.

        See: `type.SpectralDataType`

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    Returns
    -------
    out: np.ndarray [shape=(..., chroma_num, time)]
        The matrix of chroma_linear

    See Also
    --------
    chroma_cqt
    chroma_octave

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract chroma_linear data

    >>> chroma_arr = af.chroma_linear(audio_arr, samplate=sr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, chroma_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(chroma_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 x_axis='time', y_axis='chroma',
    >>>                 title='Chroma_linear')
    >>> fig.colorbar(img, ax=ax)
    """
    spec_obj = Spectrogram(num=chroma_num, radix2_exp=radix2_exp, samplate=samplate,
                           low_fre=low_fre, high_fre=high_fre, window_type=window_type,
                           slide_length=slide_length,
                           filter_bank_type=SpectralFilterBankType.CHROMA,
                           style_type=style_type, normal_type=normal_type,
                           data_type=data_type)

    spec_arr = spec_obj.spectrogram(X)
    return spec_arr


def chroma_octave(X, chroma_num=12, radix2_exp=12, samplate=32000, low_fre=note_to_hz('C1'), high_fre=16000.,
                  window_type=WindowType.HANN, slide_length=1024, data_type=SpectralDataType.POWER,
                  style_type=SpectralFilterBankStyleType.SLANEY, normal_type=SpectralFilterBankNormalType.NONE):
    """
    Octave chromagram

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    chroma_num: int
        Number of chroma bins to generate.

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    high_fre: float or None
        Highest frequency. Default is `16000(samplate/2)`.

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    data_type: SpectralDataType
        Spectrogram data type.

        See: `type.SpectralDataType`

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

        see: `type.SpectralFilterBankStyleType`

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    Returns
    -------
    out: np.ndarray [shape=(..., chroma_num, time)]
        The matrix of chroma_octave

    See Also
    --------
    chroma_linear
    chroma_cqt

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract chroma_octave data

    >>> chroma_arr = af.chroma_octave(audio_arr, samplate=sr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, chroma_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(chroma_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 x_axis='time', y_axis='chroma',
    >>>                 title='Chroma_octave')
    >>> fig.colorbar(img, ax=ax)
    """
    spec_obj = Spectrogram(num=chroma_num, radix2_exp=radix2_exp, samplate=samplate,
                           low_fre=low_fre, high_fre=high_fre, window_type=window_type,
                           slide_length=slide_length, data_type=data_type,
                           filter_bank_type=SpectralFilterBankType.OCTAVE_CHROMA,
                           style_type=style_type, normal_type=normal_type)

    spec_arr = spec_obj.spectrogram(X)
    return spec_arr


def chroma_cqt(X, chroma_num=12, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
               factor=1., thresh=0.01,
               window_type=WindowType.HANN, slide_length=None,
               normal_type=SpectralFilterBankNormalType.AREA,
               is_scale=True):
    """
    Constant-Q chromagram

    .. Note:: We recommend using the :ref:`CQT <transforms/cqt:CQT>` class, you can use it more flexibly and efficiently.

    Parameters
    ----------
    X: np.ndarray [shape=(..., n)]
        audio time series.

    chroma_num: int
        Number of chroma bins to generate.

    num: int
        Number of frequency bins to generate, starting at `low_fre`.

        Usually: ``num = octave * bin_per_octave``, `default: 84 (7 * 12)`

    samplate: int:
        Sampling rate of the incoming audio.

    low_fre: float
        Lowest frequency. `default: 32.703(C1)`

    bin_per_octave: int
        Number of bins per octave.

    factor: float
        Factor value

    thresh: float
        Thresh value

    window_type: WindowType
        Window type for each frame.

        See: `type.WindowType`

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

        See: `type.SpectralFilterBankNormalType`

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(..., chroma_num, time)]
        The matrix of CHROMA_CQT

    See Also
    --------
    chroma_linear
    chroma_octave

    Examples
    --------

    Read 220Hz audio data

    >>> import audioflux as af
    >>> audio_path = af.utils.sample_path('220')
    >>> audio_arr, sr = af.read(audio_path)

    Extract chroma_cqt data

    >>> chroma_arr = af.chroma_cqt(audio_arr, samplate=sr)

    Show spectrogram plot

    >>> import matplotlib.pyplot as plt
    >>> from audioflux.display import fill_spec
    >>> import numpy as np
    >>>
    >>> # calculate x-coords
    >>> audio_len = audio_arr.shape[-1]
    >>> x_coords = np.linspace(0, audio_len/sr, chroma_arr.shape[-1] + 1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> img = fill_spec(chroma_arr, axes=ax,
    >>>                 x_coords=x_coords,
    >>>                 x_axis='time', y_axis='chroma',
    >>>                 title='Chroma_cqt')
    >>> fig.colorbar(img, ax=ax)
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=0, thresh=thresh)
    cqt_arr = cqt_obj.cqt(X)
    power_arr = cqt_arr ** 2
    chroma_arr = cqt_obj.chroma(power_arr, chroma_num=chroma_num, data_type=SpectralDataType.POWER)
    return chroma_arr
