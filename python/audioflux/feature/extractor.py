import warnings
from typing import Iterable
from collections import defaultdict

import numpy as np
from audioflux import BFT, NSGT, CWT, PWT, CQT, ST, FST, DWT, WPT
from audioflux.feature.xxcc import XXCC
from audioflux.feature.spectral import Spectral
from audioflux.feature.deconv import Deconv
from audioflux.type import SpectralFilterBankScaleType, CepstralRectifyType, WindowType, \
    SpectralFilterBankStyleType, SpectralFilterBankNormalType, SpectralDataType, \
    NSGTFilterBankType, WaveletContinueType, WaveletDiscreteType
from audioflux.utils import note_to_hz, check_audio

__all__ = ['FeatureExtractor']


class FeatureResult(object):

    def __init__(self, name):
        self.name = name
        self._r = defaultdict(list)

    def items(self):
        return self._r.items()

    def __len__(self):
        return len(self._r)

    def __getitem__(self, item):
        return self._r[item]

    def __repr__(self):
        return f'<{self.name}: {self._r.__repr__()}>'

    def __str__(self):
        return self.__repr__()


class FeatureExtractor(object):
    """
    Batch feature extraction

    Parameters
    ----------
    transforms: list or str
        Transform type.

        Supported: bft/nsgt/cwt/pwt/cqt/st/fst/dwt/wpt

    num: int or None
        Number of frequency bins to generate.

        * only `BFT`/`NSGT`/`CWT`/`PWT` available
        * cqt is a fixed value of 84
        * dwt see: `DWT` default value
        * wpt see: `WPT` default value

    radix2_exp: int
        ``fft_length=2**radix2_exp``

    samplate: int
        Sampling rate of the incoming audio.

    low_fre: float or None
        Lowest frequency.

        * only `BFT`/`NSGT`/`CWT`/`PWT` available

    high_fre: float or None
        Highest frequency.

        * only `BFT`/`NSGT`/`CWT`/`PWT` available

    bin_per_octave: int
        Number of bins per octave.

        * only the Octave type of `BFT`/`NSGT`/`CWT`/`PWT` and `CQT` is available

    slide_length: int or None
        Window sliding length.

        * only `BFT`/`CQT` available

    scale_type: SpectralFilterBankScaleType
        Spectral filter bank type. It determines the type of spectrogram.

        * only `BFT`/`NSGT`/`CWT`/`PWT` available

        See: `type.SpectralFilterBankScaleType`

    wavelet_type: WaveletContinueType
        Wavelet Type

        * only `CWT` available

        See: `type.WaveletContinueType`

    Examples
    --------

    Get a 880Hz's audio file

    >>> import audioflux as af
    >>> sample_path = af.utils.sample_path('880')
    >>> audio_arr, sr = af.read(sample_path)

    Create FeatureExtractor object and extract spectrogram

    >>> from audioflux.type import SpectralFilterBankScaleType
    >>> fa_obj = af.FeatureExtractor(transforms=['bft', 'cwt', 'cqt'], samplate=sr, radix2_exp=12,
    >>>                              scale_type=SpectralFilterBankScaleType.OCTAVE)
    >>> spec_result = fa_obj.spectrogram(audio_arr, is_continue=True)

    Extract spectral/xxcc/deconv

    >>> spectral_result = fa_obj.spectral(spec_result, spectral='flux',
    >>>                                   spectral_kw={'is_positive': True})
    >>> xxcc_result = fa_obj.xxcc(spec_result, cc_num=13)
    >>> deconv_result = fa_obj.deconv(spec_result)
    """
    _T_BFT = 'bft'
    _T_NSGT = 'nsgt'
    _T_CWT = 'cwt'
    _T_PWT = 'pwt'
    _T_CQT = 'cqt'
    _T_ST = 'st'
    _T_FST = 'fst'
    _T_DWT = 'dwt'
    _T_WPT = 'wpt'

    def __init__(self, transforms, num=None, radix2_exp=12, samplate=32000,
                 low_fre=None, high_fre=None, bin_per_octave=12, slide_length=None,
                 scale_type=SpectralFilterBankScaleType.LINEAR,
                 wavelet_type=WaveletContinueType.MORSE):

        if radix2_exp > 23:
            raise ValueError(f'radix2_exp={radix2_exp} must be less than or equal to 23')

        if num is None:
            if scale_type in (SpectralFilterBankScaleType.LINEAR,
                              SpectralFilterBankScaleType.LINSPACE):
                num = (1 << radix2_exp) // 2 + 1
            elif scale_type in (SpectralFilterBankScaleType.MEL,
                                SpectralFilterBankScaleType.BARK,
                                SpectralFilterBankScaleType.ERB):
                num = 128
            elif scale_type in (SpectralFilterBankScaleType.OCTAVE,
                                SpectralFilterBankScaleType.LOG):
                num = 84

        self.fft_length = 1 << radix2_exp
        self.num = num
        self.radix2_exp = radix2_exp
        self.samplate = samplate
        self.low_fre = low_fre
        self.high_fre = high_fre
        self.bin_per_octave = bin_per_octave
        self.slide_length = slide_length
        self.scale_type = scale_type
        self.wavelet_type = wavelet_type

        self._transforms = {}
        if isinstance(transforms, Iterable):
            for name in transforms:
                _obj = self._create_obj(name)
                self._transforms[name] = _obj
        elif isinstance(transforms, str):
            _obj = self._create_obj(transforms)
            self._transforms[transforms] = _obj
        else:
            raise ValueError(f'transforms must be type of list or str')

    def _create_obj(self, name):
        name = name.lower()
        if name == self._T_BFT:
            return BFT(num=self.num, radix2_exp=self.radix2_exp, samplate=self.samplate,
                       low_fre=self.low_fre, high_fre=self.high_fre,
                       bin_per_octave=self.bin_per_octave,
                       window_type=WindowType.HANN, slide_length=self.slide_length,
                       scale_type=self.scale_type,
                       style_type=SpectralFilterBankStyleType.SLANEY,
                       normal_type=SpectralFilterBankNormalType.NONE,
                       data_type=SpectralDataType.MAG,
                       is_reassign=False, is_temporal=False)
        elif name == self._T_NSGT:
            return NSGT(num=self.num, radix2_exp=self.radix2_exp, samplate=self.samplate,
                        low_fre=self.low_fre, high_fre=self.high_fre,
                        bin_per_octave=self.bin_per_octave, min_len=3,
                        nsgt_filter_bank_type=NSGTFilterBankType.EFFICIENT,
                        scale_type=self.scale_type,
                        style_type=SpectralFilterBankStyleType.SLANEY,
                        normal_type=SpectralFilterBankNormalType.BAND_WIDTH)
        elif name == self._T_CWT:
            return CWT(num=self.num, radix2_exp=self.radix2_exp, samplate=self.samplate,
                       low_fre=self.low_fre, high_fre=self.high_fre,
                       bin_per_octave=self.bin_per_octave,
                       wavelet_type=self.wavelet_type,
                       scale_type=self.scale_type,
                       is_padding=True)
        elif name == self._T_PWT:
            return PWT(num=self.num, radix2_exp=self.radix2_exp, samplate=self.samplate,
                       low_fre=self.low_fre, high_fre=self.high_fre,
                       bin_per_octave=self.bin_per_octave,
                       scale_type=self.scale_type,
                       style_type=SpectralFilterBankStyleType.SLANEY,
                       normal_type=SpectralFilterBankNormalType.NONE,
                       is_padding=True)
        elif name == self._T_CQT:
            return CQT(num=84, samplate=self.samplate, low_fre=note_to_hz('C1'),
                       bin_per_octave=self.bin_per_octave,
                       window_type=WindowType.HANN, slide_length=self.slide_length,
                       normal_type=SpectralFilterBankNormalType.AREA)
        elif name == self._T_ST:
            # index: 1-2047
            return ST(radix2_exp=self.radix2_exp, samplate=self.samplate)
        elif name == self._T_FST:
            return FST(radix2_exp=self.radix2_exp, samplate=self.samplate)
        elif name == self._T_DWT:
            return DWT(num=None, radix2_exp=self.radix2_exp, samplate=self.samplate,
                       wavelet_type=WaveletDiscreteType.SYM, t1=4, t2=0)
        elif name == self._T_WPT:
            return WPT(num=None, radix2_exp=self.radix2_exp, samplate=self.samplate,
                       wavelet_type=WaveletDiscreteType.SYM, t1=4, t2=0)
        else:
            raise ValueError(f'transform name={name} is not supported')

    def _continue_spec(self, spec_fn, data_arr):
        data_len = data_arr.shape[-1]
        cur_spec_arr = None
        per_time_len = 0
        ret_arr = []
        for i, sample_idx in enumerate(range(0, data_len, self.fft_length // 2)):
            sample_arr = data_arr[..., sample_idx:sample_idx + self.fft_length].copy()
            if sample_arr.shape[-1] != self.fft_length:
                break
            cur_spec_arr = spec_fn(sample_arr)

            n_time = cur_spec_arr.shape[-1]
            if n_time < 4:
                raise ValueError(f'The length={n_time} of the time dimension of the spectrogram must be greater than 4')
            per_time_len = n_time // 4

            start_idx = 0 if i == 0 else per_time_len

            valid_spec_arr = cur_spec_arr[..., :, start_idx:per_time_len * 3]
            ret_arr.append(valid_spec_arr)
        ret_arr = np.concatenate(ret_arr, axis=-1)

        # The last part of the last package of stitching
        if cur_spec_arr is not None and per_time_len:
            ret_arr = np.concatenate((ret_arr, cur_spec_arr[..., :, per_time_len * 3:]), axis=-1)
        return ret_arr

    def _check_name(self, name):
        if name not in (self._T_BFT, self._T_NSGT, self._T_CWT, self._T_PWT,
                        self._T_CQT, self._T_ST, self._T_FST, self._T_DWT, self._T_WPT):
            raise ValueError(f'name={name} is not supported')
        return True

    def spectrogram(self, data_arr, is_continue=False):
        """
        Get the spectrogram of transforms

        Parameters
        ----------
        data_arr: np.ndarray [shape=(n,)] or list(np.ndarray [shape=(n,)])
            Input audio datas

        is_continue: bool
            Calculate continuous 2D(nsgt/cwt/pwt/st/fst/dwt/wpt) spectrogram.

            Calculate the cwt every fft/2, then divide the x-axis of each cwt into four parts, and
            take the middle two parts for splicing (the first part will be spliced with the first
            part of the first cwt, and the tail will be spliced with the last part of the last cwt).

            If True, then calculate continuous 2D spectrogram.
            If False, the 2D spectrogram can only calculate fft_length(2**radix2_exp) data

        Returns
        -------
        out: FeatureResult
            spectrogram data
        """
        data_arr = np.asarray(data_arr, dtype=np.float32, order='C')
        check_audio(data_arr, is_mono=False)

        result = FeatureResult('spectrogram')
        for name, _obj in self._transforms.items():
            if name == self._T_BFT:
                spec_fn = _obj.bft
            elif name == self._T_NSGT:
                spec_fn = _obj.nsgt
            elif name == self._T_CWT:
                spec_fn = _obj.cwt
            elif name == self._T_PWT:
                spec_fn = _obj.pwt
            elif name == self._T_CQT:
                spec_fn = _obj.cqt
            elif name == self._T_ST:
                spec_fn = _obj.st
            elif name == self._T_FST:
                spec_fn = _obj.fst
            elif name == self._T_DWT:
                spec_fn = lambda x: _obj.dwt(x)[1]
            elif name == self._T_WPT:
                spec_fn = lambda x: _obj.wpt(x)[1]
            else:
                raise ValueError(f'name={name} is not supported')

            _need_continue = is_continue
            if _need_continue and name in (self._T_BFT, self._T_CQT):
                _need_continue = False

            if _need_continue:
                spec_arr = self._continue_spec(spec_fn, data_arr)
            else:
                spec_arr = spec_fn(data_arr)
            result[name].append(spec_arr)
        return result

    def xxcc(self, spec_result, cc_num=13, rectify_type=CepstralRectifyType.LOG,
             spec_convert=np.abs):
        """
        cepstral coefficients

        See: `XXCC`

        Parameters
        ----------
        spec_result: FeatureResult
            spectrogram result

        cc_num: int
            xxcc num, usually set to 13, 20 or 40

        rectify_type: CepstralRectifyType
            rectify type

        spec_convert: function
            Operations on spectrogram data

        Returns
        -------
        out: FeatureResult
            xxcc data
        """
        result = FeatureResult(name='xxcc')
        for name, spec_list in spec_result.items():
            if name not in self._transforms:
                raise ValueError(f'name={name} not found in transforms')

            _obj = self._transforms[name]
            xxcc_obj = XXCC(_obj.num)
            for spec in spec_list:
                n_time = spec.shape[-1]
                xxcc_obj.set_time_length(n_time)
                spec = spec_convert(spec)
                xx_arr = xxcc_obj.xxcc(spec, cc_num=cc_num,
                                       rectify_type=rectify_type)
                result[name].append(xx_arr)
        return result

    def deconv(self, spec_result, spec_convert=np.abs):
        """
        Deconv feature

        See: `Deconv`

        Parameters
        ----------
        spec_result: FeatureResult
            spectrogram result

        spec_convert: function
            Operations on spectrogram data

        Returns
        -------
        out: FeatureResult
            deconv data
        """
        result = FeatureResult('deconv')
        for name, spec_list in spec_result.items():
            if name not in self._transforms:
                raise ValueError(f'name={name} not found in transforms')
            if not spec_list:
                continue

            _obj = self._transforms[name]
            deconv_obj = Deconv(_obj.num)
            n_time = spec_list[0].shape[-1]
            deconv_obj.set_time_length(n_time)
            for spec in spec_list:
                spec = spec_convert(spec)
                deconv_arr = deconv_obj.deconv(spec)
                result[name].append(deconv_arr)
        return result

    def spectral(self, spec_result, spectral, spectral_kw=None, spec_convert=np.abs):
        """
        spectral related features

        Parameters
        ----------
        spec_result: FeatureResult
            spectrogram result

        spectral: str
            spectral feature name

            See: `Spectral`

        spectral_kw: dict or None
            spectral feature parameters

            See: `Spectral`

        spec_convert: function
            Operations on spectrogram data

        Returns
        -------
        out: FeatureResult
            spectral data
        """
        result = FeatureResult(f'spectral.{spectral}')
        for name, spec_list in spec_result.items():
            if name not in self._transforms:
                raise ValueError(f'name={name} not found in transforms')
            if not spec_list:
                continue

            _obj = self._transforms[name]
            spectral_obj = Spectral(num=_obj.num,
                                    fre_band_arr=_obj.get_fre_band_arr())
            if not hasattr(spectral_obj, spectral):
                raise ValueError(f'spectral={spectral} function is not found')
            n_time = spec_list[0].shape[-1]
            spectral_obj.set_time_length(n_time)
            for spec in spec_list:
                spec = spec_convert(spec)
                spectral_kw = {} if spectral_kw is None else spectral_kw
                spectral_arr = getattr(spectral_obj, spectral)(spec, **spectral_kw)
                result[name].append(spectral_arr)
        return result
