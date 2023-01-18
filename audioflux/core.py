import numpy as np
from audioflux import BFT, XXCC, CQT
from audioflux.spectrogram import Spectrogram
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
                       is_reassign=False, result_type=0):
    """
    Short-time Fourier transform (Linear/STFT)

    It is a Fourier-related transform used to determine
    the sinusoidal frequency and phase content of
    local sections of a signal as it changes over time.

    .. Note:: If you want to set more parameters, use the `BFT` class

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

    is_reassign: bool
        Whether to use reassign.

    result_type: int，0 or 1
        - If `0`, then the result is a matrix of complex numbers.
        - If `1`, then the result is a matrix of real numbers.

    Returns
    -------
    out: np.ndarray [shape=(fre, time)]
        The matrix of Linear(STFT)

    See Also
    --------
    mel_spectrogram
    bark_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT
    """
    if num is None:
        num = (1 << radix2_exp) / 2 + 1
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.LINEAR, style_type=style_type,
                  data_type=data_type, is_reassign=is_reassign)
    spec_arr = bft_obj.bft(X, result_type=result_type)
    return spec_arr


def mel_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                    low_fre=0., high_fre=None,
                    window_type=WindowType.HANN, slide_length=None,
                    style_type=SpectralFilterBankStyleType.SLANEY,
                    normal_type=SpectralFilterBankNormalType.NONE,
                    data_type=SpectralDataType.POWER,
                    is_reassign=False, result_type=0):
    """
    Mel-scale spectrogram.

    .. Note:: If you want to set more parameters, use the `BFT` class

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

    is_reassign: bool
        Whether to use reassign.

    result_type: int，0 or 1
        - If `0`, then the result is a matrix of complex numbers.
        - If `1`, then the result is a matrix of real numbers.

    Returns
    -------
    out: np.ndarray [shape=(fre, time)]
        The matrix of MEL

    See Also
    --------
    linear_spectrogram
    bark_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.MEL, style_type=style_type,
                  normal_type=normal_type, data_type=data_type,
                  is_reassign=is_reassign)
    spec_arr = bft_obj.bft(X, result_type=result_type)
    return spec_arr


def bark_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                     low_fre=None, high_fre=None,
                     window_type=WindowType.HANN, slide_length=None,
                     style_type=SpectralFilterBankStyleType.SLANEY,
                     normal_type=SpectralFilterBankNormalType.NONE,
                     data_type=SpectralDataType.POWER,
                     is_reassign=False, result_type=0):
    """
    Bark-scale spectrogram.

    .. Note:: If you want to set more parameters, use the `BFT` class

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

    is_reassign: bool
        Whether to use reassign.

    result_type: int，0 or 1
        - If `0`, then the result is a matrix of complex numbers.
        - If `1`, then the result is a matrix of real numbers.

    Returns
    -------
    out: np.ndarray [shape=(fre, time)]
        The matrix of BARK

    See Also
    --------
    linear_spectrogram
    mel_spectrogram
    erb_spectrogram

    BFT
    NSGT
    CWT
    PWT
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.BARK, style_type=style_type,
                  normal_type=normal_type, data_type=data_type,
                  is_reassign=is_reassign)
    spec_arr = bft_obj.bft(X, result_type=result_type)
    return spec_arr


def erb_spectrogram(X, num=128, radix2_exp=12, samplate=32000,
                    low_fre=None, high_fre=None,
                    window_type=WindowType.HANN, slide_length=None,
                    style_type=SpectralFilterBankStyleType.SLANEY,
                    normal_type=SpectralFilterBankNormalType.NONE,
                    data_type=SpectralDataType.POWER,
                    is_reassign=False, result_type=0):
    """
    Erb-scale spectrogram.

    .. Note:: If you want to set more parameters, use the `BFT` class

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int or None
        Window sliding length.

        If `slide_length` is None, then ``slide_length = fft_length / 4``

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    data_type: SpectralDataType
        Spectrogram data type.

        It cat be set to mag or power. If you needs `db` type,
        you can set `power` type and then call the `power_to_db` method.

    is_reassign: bool
        Whether to use reassign.

    result_type: int，0 or 1
        - If `0`, then the result is a matrix of complex numbers.
        - If `1`, then the result is a matrix of real numbers.

    Returns
    -------
    out: np.ndarray [shape=(fre, time)]
        The matrix of ERB

    See Also
    --------
    linear_spectrogram
    mel_spectrogram
    bark_spectrogram

    BFT
    NSGT
    CWT
    PWT
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.ERB, style_type=style_type,
                  normal_type=normal_type, data_type=data_type,
                  is_reassign=is_reassign)
    spec_arr = bft_obj.bft(X, result_type=result_type)
    return spec_arr


def mfcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Mel-frequency cepstral coefficients (MFCCs)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio time series.

    cc_num: int
        number of MFCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    num: int
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

    Returns
    -------
    out: np.ndarray [shape=(cc_num, time)]
        The matrix of MFCCs

    See Also
    --------
    bfcc
    gtcc
    cqcc

    BFT
    XXCC
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.MEL,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    return xxcc_arr


def bfcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Bark-frequency cepstral coefficients (BFCCs)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio time series.

    cc_num: int
        number of BFCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    num: int
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

    Returns
    -------
    out: np.ndarray [shape=(cc_num, time)]
        The matrix of BFCCs

    See Also
    --------
    mfcc
    gtcc
    cqcc

    BFT
    XXCC
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.BARK,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    return xxcc_arr


def gtcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         num=128, radix2_exp=12, samplate=32000, slide_length=None,
         low_fre=None, high_fre=None, window_type=WindowType.HANN):
    """
    Gammatone cepstral coefficients (GTCCs)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio time series.

    cc_num: int
        number of GTCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    num: int
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

    Returns
    -------
    out: np.ndarray [shape=(cc_num, time)]
        The matrix of GTCCs

    See Also
    --------
    mfcc
    bfcc
    cqcc

    BFT
    XXCC
    """
    bft_obj = BFT(num=num, radix2_exp=radix2_exp, samplate=samplate,
                  low_fre=low_fre, high_fre=high_fre,
                  window_type=window_type, slide_length=slide_length,
                  scale_type=SpectralFilterBankScaleType.ERB,
                  style_type=SpectralFilterBankStyleType.GAMMATONE,
                  normal_type=SpectralFilterBankNormalType.AREA,
                  data_type=SpectralDataType.POWER, is_reassign=False)
    spec_arr = bft_obj.bft(X)
    print(spec_arr.shape)

    xxcc_obj = XXCC(num=bft_obj.num)
    xxcc_obj.set_time_length(time_length=spec_arr.shape[1])
    xxcc_arr = xxcc_obj.xxcc(spec_arr, cc_num=cc_num, rectify_type=rectify_type)
    return xxcc_arr


def cqcc(X, cc_num=13, rectify_type=CepstralRectifyType.LOG,
         num=84, samplate=32000, low_fre=note_to_hz('C1'), slide_length=None,
         bin_per_octave=12, window_type=WindowType.HANN,
         normal_type=SpectralFilterBankNormalType.AREA,
         is_scale=True):
    """
    Constant-Q cepstral coefficients (CQCCs)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
        audio time series.

    cc_num: int
        number of GTCC to return.

    rectify_type: CepstralRectifyType
        cepstral rectify type

    num: int
        Number of mel frequency bins to generate, starting at `low_fre`.

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

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(cc_num, time)]
        The matrix of GTCCs

    See Also
    --------
    mfcc
    bfcc
    gtcc

    CQT
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale)
    spec_arr = cqt_obj.cqt(X)
    print(spec_arr.shape)
    power_arr = np.abs(spec_arr) ** 2
    cqcc_arr = cqt_obj.cqcc(power_arr, cc_num=cc_num, rectify_type=rectify_type)
    return cqcc_arr


def cqt(X, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
        factor=1., thresh=0.01,
        window_type=WindowType.HANN, slide_length=None,
        normal_type=SpectralFilterBankNormalType.AREA,
        is_scale=True):
    """
    Constant-Q transform (CQT)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(num, time)]
        The matrix of CQT
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=0., thresh=thresh)
    cqt_arr = cqt_obj.cqt(X)
    return cqt_arr


def vqt(X, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
        factor=1., beta=0.5, thresh=0.01,
        window_type=WindowType.HANN, slide_length=None,
        normal_type=SpectralFilterBankNormalType.AREA,
        is_scale=True):
    """
    Variable-Q transform (VQT)

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(num, time)]
        The matrix of CQT
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=beta, thresh=thresh)
    cqt_arr = cqt_obj.cqt(X)
    return cqt_arr


def chroma_linear(X, chroma_num=12, radix2_exp=12, samplate=32000, low_fre=0., high_fre=16000.,
                  window_type=WindowType.HANN, slide_length=1024,
                  style_type=SpectralFilterBankStyleType.SLANEY,
                  normal_type=SpectralFilterBankNormalType.NONE):
    """
    Linear(STFT) chromagram

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int
        Window sliding length.

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    Returns
    -------
    out: np.ndarray [shape=(chroma_num, time)]
        The matrix of chroma_linear

    See Also
    --------
    chroma_cqt
    chroma_octave
    """
    spec_obj = Spectrogram(num=chroma_num, radix2_exp=radix2_exp, samplate=samplate,
                           low_fre=low_fre, high_fre=high_fre, window_type=window_type,
                           slide_length=slide_length,
                           filter_bank_type=SpectralFilterBankType.CHROMA,
                           filter_style_type=style_type, filter_normal_type=normal_type)

    spec_arr = spec_obj.spectrogram(X)
    return spec_arr


def chroma_octave(X, chroma_num=12, radix2_exp=12, samplate=32000, low_fre=0., high_fre=16000.,
               window_type=WindowType.HANN, slide_length=1024,
               style_type=SpectralFilterBankStyleType.SLANEY, normal_type=SpectralFilterBankNormalType.NONE):
    """
    Octave chromagram

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    slide_length: int
        Window sliding length.

    style_type: SpectralFilterBankStyleType
        Spectral filter bank style type. It determines the bank type of window.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    Returns
    -------
    out: np.ndarray [shape=(chroma_num, time)]
        The matrix of chroma_octave

    See Also
    --------
    chroma_linear
    chroma_cqt
    """
    spec_obj = Spectrogram(num=chroma_num, radix2_exp=radix2_exp, samplate=samplate,
                           low_fre=low_fre, high_fre=high_fre, window_type=window_type,
                           slide_length=slide_length,
                           filter_bank_type=SpectralFilterBankType.OCTAVE_CHROMA,
                           filter_style_type=style_type, filter_normal_type=normal_type)

    spec_arr = spec_obj.spectrogram(X)
    return spec_arr


def chroma_cqt(X, chroma_num=12, num=84, samplate=32000, low_fre=note_to_hz('C1'), bin_per_octave=12,
               factor=1., beta=0., thresh=0.01,
               window_type=WindowType.HANN, slide_length=None,
               normal_type=SpectralFilterBankNormalType.AREA,
               is_scale=True):
    """
    Constant-Q chromagram

    Parameters
    ----------
    X: np.ndarray [shape=(n,)]
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

    beta: float
        Beta value

    thresh: float
        Thresh value

    window_type: WindowType
        Window type for each frame.

    slide_length: int
        Window sliding length.

    normal_type: SpectralFilterBankNormalType
        Spectral filter normal type. It determines the type of normalization.

    is_scale: bool
        Whether to use scale.

    Returns
    -------
    out: np.ndarray [shape=(chroma_num, time)]
        The matrix of CHROMA_CQT

    See Also
    --------
    chroma_linear
    chroma_octave
    """
    cqt_obj = CQT(num=num, samplate=samplate, low_fre=low_fre, slide_length=slide_length,
                  bin_per_octave=bin_per_octave, window_type=window_type,
                  normal_type=normal_type, is_scale=is_scale,
                  factor=factor, beta=beta, thresh=thresh)
    cqt_arr = cqt_obj.cqt(X)
    chroma_arr = cqt_obj.chroma(cqt_arr, chroma_num=chroma_num)
    return chroma_arr
