from enum import Enum

__all__ = [
    'WindowType',
    'FilterBandType',
    'SpectralDataType',
    'SpectralFilterBankType',
    'SpectralFilterBankScaleType',
    'SpectralFilterBankStyleType',
    'SpectralFilterBankNormalType',
    'SpectralNoveltyMethodType',
    'SpectralNoveltyDataType',
    'ChromaDataNormalType',
    'CepstralRectifyType',
    'CepstralEnergyType',
    'PaddingPositionType',
    'PaddingModeType',
    'WaveletContinueType',
    'WaveletDiscreteType',

    'get_wavelet_default_gamma_beta'
]


class WindowType(Enum):
    """
    Window Type

    Attributes
    ----------
    HANN:
        :math:`\quad w(n)=0.5\\left(  1-\cos \\left(2\pi \cfrac n{N} \\right)\\right) , 0 \le n \le N`

    HAMM:
        :math:`\quad w(n)=0.54 - 0.46\cos \\left(2\pi \cfrac n{N} \\right)\ , 0 \le n \le N`

    BLACKMAN:
        :math:`\quad w(n)=0.42 - 0.5\cos \\left(2\pi \cfrac n{N-1} \\right)\ +0.08\cos\\left(  4\pi\cfrac n{N-1} \\right)  , 0 \le n \le M-1`

        M=N/2,(N+1/2) 当N为偶数、奇数时

    KAISER:
        :math:`\quad w(n)=\cfrac {I_0 \\left(  \beta \sqrt{1- \\left( { \cfrac {n-N/2}{N/2} }  \\right)^2  }  \\right)  } {I_0(\beta)}  , 0 \le n \le N`

        | :math:`I_0(\\beta)` 为零阶第一类修正贝塞尔函数，可有下面公式级数计算，
        | :math:`I_0(\beta)=1+\sum_{k=1}^{\infty} \\left[ \cfrac1{k!} \\left (\cfrac  \beta 2 \\right)^k \\right] ^2`，一般取15项左右
        | :math:`\\beta=\\begin {cases} 0.1102(\\alpha_s-8.7),  & \\alpha_s>50 \\\\ 0.5842(\\alpha_s-21)^{0.4}+0.07886(\\alpha_s-21), & 50 \ge \\alpha_s \ge 21 \\\\   0, & \\alpha_s <21 \end {cases} \quad \\alpha_s` 为旁瓣衰减dB
        | 默认 :math:`\\beta=5`

    BARTLETT:
        :math:`\quad w(n)=\\begin{cases}  \cfrac{2n}N, & 0 \le n \le \cfrac N{2}  \\  2-\cfrac{2n}N， & \cfrac N{2} \le n \le N  \end{cases}`

        Bartlett和Triang非常相似，Bartlett首尾处为0，Triang不为0

    TRIANG:
        | :math:`N` 为奇数时，
        | :math:`\quad w(n)=\\begin{cases}  \cfrac{2n}{N+1}, & 1 \le n \le \cfrac {N+1}{2}  \\\\  2-\cfrac{2n}{N+1}， & \cfrac {N+1}2 \le n \le N  \end{cases}`
        | :math:`N` 为偶数时，
        | :math:`\quad w(n)=\\begin{cases}  \cfrac{2n-1}N, & 1 \le n \le \cfrac N{2}  \\\\  2-\cfrac{2n-1}N， & \cfrac{N}2 +1 \le n \le N  \end{cases}`
    FLATTOP:
        :math:`\quad w(n)=a_0-a_1\cos\\left( \cfrac{2\pi n}{N-1} \\right) + a_2\cos\\left( \cfrac{4\pi n}{N-1} \\right) -a_3\cos\\left( \cfrac{6\pi n}{N-1} \\right) +a_4\cos\\left( \cfrac{8\pi n}{N-1} \\right)`

        :math:`\\begin{cases}  a_0 =0.21557895 \\\\ a_1= 0.41663158 \\\\ a_2= 0.277263158 \\\\ a_3= 0.083578947 \\\\ a_4= 0.006947368  \end{cases}`

    GAUSS:
        :math:`\quad w(n)=e^{-n^2/2\sigma^2} =e^{ -\\frac12 \\left(  \\alpha \\frac{n}{ (N-1)/2 }  \\right)^2 }  , \qquad -(N-1)/2 \le n \le (N-1)/2`

        x使用linspace产生，使用窗长度N，窗首尾强制设为0

    TUKEY:
        :math:`w(x)= \\begin{cases} \\frac1{2}(1+ \cos(\\frac{2\pi}{\\alpha} [x-\\alpha/2] ) ) ,& 0 \le x \le \cfrac\\alpha{2} \\\\ 1,       & \cfrac \\alpha{2} \le x <1-\cfrac\\alpha{2}   \\\\  \\frac1{2}(1+ \cos(\\frac{2\pi}{\\alpha} [x-1+\\alpha/2] ) ) ,& 1-\cfrac\\alpha{2}  \le x \le 1   \end{cases}`

        | x使用linspace产生，使用窗长度N
        | 默认 :math:`\\alpha=0.5`
    """
    RECT = 0
    HANN = 1
    HAMM = 2

    BLACKMAN = 3
    KAISER = 4

    BARTLETT = 5
    TRIANG = 6

    FLATTOP = 7
    GAUSS = 8

    BLACKMAN_HARRIS = 9
    BLACKMAN_NUTTALL = 10
    BARTLETT_HANN = 11

    BOHMAN = 12

    TUKEY = 13  # tapered cosine


class FilterBandType(Enum):
    """
    Filter Band Type
    """
    LOW_PASS = 0
    HIGH_PASS = 1

    BAND_PASS = 2
    BAND_STOP = 3  # Rejection

    # ALL_PASS = 4


class SpectralDataType(Enum):
    """
    Spectral Data Type
    """
    #: power type
    POWER = 0
    #: magnitude type
    MAG = 1


class SpectralFilterBankScaleType(Enum):
    """
    Spectral Filter Bank Scale Type
    """
    #: Short-time Fourier transform spectrogram.
    LINEAR = 0
    #: Linspace-scale spectrogram.
    LINSPACE = 1

    #: Mel-scale spectrogram.
    MEL = 2
    #: Bark-scale spectrogram.
    BARK = 3
    #: Erb-scale spectrogram.
    ERB = 4

    #: Octave-scale spectrogram.
    OCTAVE = 5  # similar Constant-Q
    #: Logarithmic-scale spectrogram.
    LOG = 6


class SpectralFilterBankType(Enum):
    """
    Spectral Filter Bank Type
    """
    #: Short-time Fourier transform spectrogram.
    LINEAR = 0
    #: Linspace-scale spectrogram.
    LINSPACE = 1

    #: Mel-scale spectrogram.
    MEL = 2
    #: Bark-scale spectrogram.
    BARK = 3
    #: Erb-scale spectrogram.
    ERB = 4

    #: Octave-scale spectrogram.
    OCTAVE = 5  # similar Constant-Q
    #: Logarithmic-scale spectrogram.
    LOG = 6

    #: Deep-scale spectrogram.
    DEEP = 7

    #: Linear(STFT) chroma.
    CHROMA = 8

    #: Octave-scale chroma. Similar cqt-chroma.
    OCTAVE_CHROMA = 9
    #: Deep-scale chroma. similar cqt-chroma.
    DEEP_CHROMA = 10


class SpectralFilterBankStyleType(Enum):
    """
    Spectral Filter Bank Style Type

    See: `type.WindowType`
    """
    SLANEY = 0  # Triang
    ETSI = 1  # Bartlett
    GAMMATONE = 2  # gammatone

    POINT = 3
    RECT = 4

    HANN = 5
    HAMM = 6

    BLACKMAN = 7
    BOHMAN = 8

    KAISER = 9
    GAUSS = 10


class SpectralFilterBankNormalType(Enum):
    """
    Spectral Filter Bank Normal Type
    """
    NONE = 0  # same hight

    AREA = 1  # normal(same hight) / same area
    BAND_WIDTH = 2


class SpectralNoveltyMethodType(Enum):
    """
    Spectral Novelty Method Type

    Attributes
    ----------
    SUB:
        :math:`sub_k(t)= s_k(t)-s_k(t-1)`
    ENTROY:
        :math:`entropy_k(t)= \log \\left( \\frac {s_k(t)}{s_k(t-1} \\right)`
    KL:
        :math:`kl_k(t)= s_k(t) \log \\left( \\frac {s_k(t)}{s_k(t-1} \\right)`
    IS:
        :math:`is_k(t)= \\frac {s_k(t)}{s_k(t-1)} - \log \\left( \\frac {s_k(t)}{s_k(t-1} \\right)-1`
    """
    SUB = 0

    ENTROY = 1
    KL = 2
    IS = 3


class SpectralNoveltyDataType(Enum):
    """
    Spectral Novelty Data Type

    * :math:`f_k=sub_k,entropy_k,\cdots,is_k \quad g_k=\log(1+\gamma f_k)`
    * :math:`v_k=f_k,g_k`

    Attributes
    ----------
    VALUE:
        :math:`\mathcal{V}(t)=\sum_{k=b_1}^{b_2}v_k(t)`,
        calculated when :math:`v_k(t) \ge \\alpha` is satisfied,
        generally :math:`\\alpha \ge 0`
    NUMBER:
        :math:`\mathcal{V}(t) =i[v_{k_{\in [b_1,b_2]}} (t) ]`,
        number statistics when :math:`v_k(t) \ge \\alpha` is satisfied,
        generally :math:`\\alpha \ge 0`

    """
    VALUE = 0
    NUMBER = 1


class ChromaDataNormalType(Enum):
    """
    Chroma Data Normal Type
    """
    NONE = 0

    MAX = 1
    MIN = 2

    P2 = 3
    P1 = 4


class CepstralRectifyType(Enum):
    """
    Cepstral Rectify Type
    """
    LOG = 0
    CUBIC_ROOT = 1


class CepstralEnergyType(Enum):
    """
    Cepstral Energy Type
    """
    REPLACE = 0
    APPEND = 1
    IGNORE = 2


class PaddingPositionType(Enum):
    """
    Padding Position Type
    """
    CENTER = 0
    RIGHT = 1
    LEFT = 2


class PaddingModeType(Enum):
    """
    Padding Mode Type
    """
    CONSTANT = 0
    REFLECT = 1
    WRAP = 2  # repeat


class WaveletContinueType(Enum):
    """
    Wavelet Continue Type

    MORSE:
        :math:`\Psi_{\beta, \gamma}(\omega) = U(\omega)a_{\\beta,\gamma}  \omega^{\\beta}e^{-\omega^{\gamma}}`

        * :math:`U(\omega)` is the step function,
          :math:`\omega` is the positive and negative sign,
          :math:`\omega^{\beta}e^{-\omega^{\gamma}}` is :math:`e^{\beta\ln\omega -\omega^{\gamma} }`
        * :math:`a_{\beta,\gamma}=\omega_0^{-\beta} e^{\omega_0^\gamma}=e^{-\beta\ln\omega_0+\omega_0^\gamma}`,
          :math:`\omega_0` represents the peak frequency of the center of the wavelet window,
          :math:`\omega_0=\\left( \\cfrac\beta\gamma \\right)^{\\frac1\gamma}` is :math:`e^{\\frac1\gamma(\ln\\beta-\ln\gamma)}` calculation
        * default: :math:`\gamma=3,\beta=20`

    MORLET:
        :math:`\Psi(\omega)=\pi^{-1/4}e^{-(\omega-\omega_0)^2/\sigma}`

        * default: :math:`\omega_0=6,\sigma=2`

    BUMP:
        :math:`\Psi(\omega)=e^{(1-\\frac{1}{1-(\omega-\omega_0)^2/\sigma^2}) }I(\omega)`

        * :math:`I(\omega)` represents the :math:`\omega` interval range mark, :math:`\|  (w-\omega_0)/\sigma  \|  \le 1`
        * default: :math:`\omega_0=5,\sigma=0.6`

    PAUL:
        :math:`\Psi(\omega)=U(\omega)\\frac{2^m}{\sqrt{m(2m-1)!}}\omega^{m}e^{-\omega}`

        * default: :math:`m=4`

    Default gamma/beta values for different wavelet_types:

    * morse: gamma=3 beta=20
    * morlet: gamma=6 beta=2
    * bump: gamma=5 beta=0.6
    * paul: gamma 4
    * dog: gamma 2 beta 2; must even
    * mexican: beta 2
    * hermit: gamma 5 beta 2
    * ricker: gamma 4
    """
    MORSE = 0
    MORLET = 1
    BUMP = 2

    PAUL = 3
    DOG = 4  # DOG
    MEXICAN = 5  # DOG order = 2

    HERMIT = 6
    RICKER = 7


class WaveletDiscreteType(Enum):
    """
    Wavelet Discrete Type

    t1/t2 settings for wavelet_type:

    * DB: t1
        * 2~10/20/30/40

    * SYM: t1
        * 2~10/20/30

    * COIF: t1
        * 1/2/3/4/5

    * FK: t1
        * 4/6/8/14/18/22

    * BIOR/DMEY: t1.t2
        * 1.1/1.3/1.5
        * 2.2/2.4/2.6/2.8
        * 3.1/3.3/3.5/3.7/3.9
        * 4.4/5.5/6.8
    """
    HAAR = 0
    DB = 1  # 2~10/20/30/40
    SYM = 2  # 2~10/20/30
    COIF = 3  # 1~5

    FK = 4  # 4/6/8/14/18/22

    # 1.1/1.3/1.5
    # 2.2/2.4/2.6/2.8
    # 3.1/3.3/3.5/3.7/3.9
    # 4.4/5.5/6.8
    BIOR = 5
    DMEY = 6


def get_wavelet_default_gamma_beta(wavelet_type):
    """
    Get default gamma/beta of wavelet

    Parameters
    ----------
    wavelet_type: WaveletContinueType
        wavelet type

        See: `WaveletContinueType`

    Notes
    -----
    morse: gamma=3 beta=20
    morlet: gamma=6 beta=2
    bump: gamma=5 beta=0.6
    paul: gamma 4
    dog: gamma 2 beta 2; must even
    mexican: beta 2
    hermit: gamma 5 beta 2
    ricker: gamma 4

    Returns
    -------
    gamma: int
    beta: int
    """
    gamma = 0
    beta = 0
    if wavelet_type == WaveletContinueType.MORSE:
        gamma = 3
        beta = 20
    elif wavelet_type == WaveletContinueType.MORLET:
        gamma = 6
        beta = 2
    elif wavelet_type == WaveletContinueType.BUMP:
        gamma = 5
        beta = 0.6
    elif wavelet_type == WaveletContinueType.PAUL:
        gamma = 4
    elif wavelet_type == WaveletContinueType.DOG:
        gamma = 2
        beta = 2
    elif wavelet_type == WaveletContinueType.MEXICAN:
        beta = 2
    elif wavelet_type == WaveletContinueType.HERMIT:
        gamma = 5
        beta = 2
    elif wavelet_type == WaveletContinueType.RICKER:
        gamma = 4
    return gamma, beta
