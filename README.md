

<img src='./image/logo.png'  width="400"  style="max-width: 100%;" > 


# audioFlux

<!--<p align="center">
<img src='./image/logo.png'  width="380"  style="max-width: 100%;" > 
</p>-->


<!-- 
[![Pypi Downloads](https://img.shields.io/pypi/dm/aubio.svg?label=Pypi%20downloads)](https://pypi.org/project/aubio/)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/aubio.svg?label=Conda%20downloads)](https://anaconda.org/conda-forge/aubio)
[![Documentation](https://readthedocs.org/projects/aubio/badge/?version=latest)](http://aubio.readthedocs.io/en/latest/?badge=latest "Latest documentation") -->
<!--![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/libAudioFlux/audioFlux)-->

<!-- shields.io -->
![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/libAudioFlux/audioFlux/build.yml?branch=master) 
![GitHub](https://img.shields.io/github/license/libAudioFlux/audioFlux)
![example branch parameter](https://github.com/libAudioFlux/audioFlux/actions/workflows/build.yml/badge.svg?branch=master)
![language](https://img.shields.io/badge/language-python%20|%20c%20-blue.svg)
[![PyPI - Version](https://img.shields.io/pypi/v/audioflux)](https://pypi.org/project/audioflux/)
[![PyPI - Python Version](https://img.shields.io/badge/python-%3E%3D3.6-brightgreen)](https://pypi.org/project/audioflux/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/audioflux.svg?label=Pypi%20downloads)](https://pypi.org/project/audioflux/)

[![Docs](https://img.shields.io/badge/Docs-passing-brightgreen)](https://audioflux.top/index.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7548288.svg)](https://doi.org/10.5281/zenodo.7548288)

<!--[![codebeat badge](https://codebeat.co/badges/0e21a344-0928-4aee-8262-be9a41fa488b)](https://codebeat.co/projects/github-com-libaudioflux-audioflux-master)
![](https://img.shields.io/badge/pod-v0.1.1-blue.svg)-->


A library for audio and music analysis, feature extraction. 


#  Table of Contents

- [Overview](#overview)
  - [Description](#description)
  - [Functionality](#functionality)
    - [Transform](#1-transform)
    - [Feature](#2-feature)
    - [MIR](#3-mir)
- [Quickstart](#quickstart)
	-  [Mel & MFCC](#mel--mfcc)
	-  [CWT & Synchrosqueezing](#cwt--synchrosqueezing)
	-  [Other examples](#other-examples)
- [Installation](#installation)
    - [Python Package Intsall](#python-package-intsall)
    - [iOS build](#ios-build)
    - [Android build](#android-build)
    - [Building from source](#building-from-source)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Citing](#citing)
- [License](#license)


## Overview 

###  Description 

In audio domain, feature extraction is particularly important for Audio Classification, Speech enhancement, Audio/Music Separation,music-information-retrieval(MIR), ASR and other audio task.

In the above tasks, **mel spectrogram** and **mfcc** features are commonly used in traditional machine-learning based on statistics and deep-learning based on neural network.

**`audioFlux`** provides systematic, comprehensive and multi-dimensional feature extraction and combination, and combines various deep learning network models to conduct research and development learning in different fields.

Can be used for deep learning, pattern recognition, signal processing, bioinformatics, statistics, finance, etc.

###  Functionality 

**`audioFlux`** is based on the design of data flow. It decouples each algorithm module structurally, and it is convenient, fast and efficient to extract features from large batches.The following are the main feature architecture diagrams, specific and detailed description view the documentation.  

<img src='./image/feature_all.png'>
<!--<img src='./feature_all.pdf'>-->

The main functions of **`audioFlux`** include **transform**, **feature** and **mir** modules. 
#### 1. Transform
In the time–frequency representation, main transform algorithm: 

- **`BFT`**&nbsp;&nbsp; -  &nbsp;&nbsp;Based Fourier Transform, similar short-time Fourier transform.
- **`NSGT`** - &nbsp; Non-Stationary Gabor Transform.
- **`CWT`**&nbsp;&nbsp; - &nbsp;&nbsp;Continuous Wavelet Transform.
- **`PWT`**&nbsp;&nbsp; - &nbsp;&nbsp;Pseudo Wavelet Transform.

<!-- &emsp -->

The above transform supports all the following frequency scale types:  

- Linear -  Short-time Fourier transform spectrogram.
- Linspace - Linspace-scale spectrogram.
- Mel - Mel-scale spectrogram.
- Bark - Bark-scale spectrogram.
- Erb - Erb-scale spectrogram.
- Octave - Octave-scale spectrogram.
- Log - Logarithmic-scale spectrogram.

The following transform are not supports multiple frequency scale types, only used as independent transform: 

- **`CQT`** -  &nbsp;&nbsp;Constant-Q Transform.
- **`VQT`** -  &nbsp;&nbsp;Variable-Q Transform. 
- **`ST`**&nbsp;&nbsp; - &nbsp;&nbsp;S-Transform/Stockwell Transform.
- **`FST`** - &nbsp;&nbsp;Fast S-Transform.
- **`DWT`** - &nbsp;&nbsp;Discrete Wavelet Transform.
- **`WPT`** - &nbsp;&nbsp;Wave Packet Transform.
- **`SWT`** - &nbsp;&nbsp;Stationary Wavelet Transform. 

Detailed transform function, description, and use view the documentation.

The *_synchrosqueezing_* or *_reassignment_* is a technique for sharpening a time-frequency representation, contains the following algorithms:

- **`reassign`** - reassign transform for `STFT`.
- **`synsq`** - reassign data use `CWT` data. 
- **`wsst`** - reassign transform for `CWT`.

#### 2. Feature
The feature module contains the following algorithms:

- **`spectral`** - Spectrum feature, supports all spectrum types.
- **`xxcc`** -  Cepstrum coefficients, supports all spectrum types.
- **`deconv`** - Deconvolution for spectrum, supports all spectrum types.
- **`chroma`** - Chroma feature, only supports `CQT` spectrum, Linear/Octave spectrum based on `BFT`. 

<!-- harmonic pitch class profiles(HPCP) -->

#### 3. MIR 
The mir module contains the following algorithms:

- **`pitch`** - YIN, STFT, etc algorithm.
- **`onset`** - Spectrum flux, novelty, etc algorithm. 
- **`hpss`** - Median filtering, NMF algorithm.


## Quickstart

### Mel & MFCC

Mel spectrogram and Mel-frequency cepstral coefficients 

```python
# Feature extraction example
import numpy as np
import audioflux as af
import matplotlib.pyplot as plt
from audioflux.display import fill_spec
from audioflux.type import SpectralFilterBankScaleType

# Get a 220Hz's audio file path
sample_path = af.utils.sample_path('220')

# Read audio data and sample rate
audio_arr, sr = af.read(sample_path)

# Extract mel spectrogram
bft_obj = af.BFT(num=128, radix2_exp=12, samplate=sr,
                 scale_type=SpectralFilterBankScaleType.MEL)
spec_arr = bft_obj.bft(audio_arr)
spec_arr = np.abs(spec_arr)

# Create XXCC object and extract mfcc
xxcc_obj = af.XXCC(bft_obj.num)
xxcc_obj.set_time_length(time_length=spec_arr.shape[1])
mfcc_arr = xxcc_obj.xxcc(spec_arr)

audio_len = audio_arr.shape[0]
fig, ax = plt.subplots()
img = fill_spec(spec_arr, axes=ax,
          x_coords=bft_obj.x_coords(audio_len),
          y_coords=bft_obj.y_coords(),
          x_axis='time', y_axis='log',
          title='Mel Spectrogram')
fig.colorbar(img, ax=ax)

fig, ax = plt.subplots()
img = fill_spec(mfcc_arr, axes=ax,
          x_coords=bft_obj.x_coords(audio_len), x_axis='time',
          title='MFCC')
fig.colorbar(img, ax=ax)

plt.show()
```

<img src='image/demo_mel.png'  width="415"  /><img src='image/demo_mfcc.png'  width="415"  />

### CWT & Synchrosqueezing

Continuous Wavelet Transform spectrogram and its corresponding synchrosqueezing reassignment spectrogram

```python
# Feature extraction example
import numpy as np
import audioflux as af
import matplotlib.pyplot as plt
from audioflux.display import fill_spec
from audioflux.type import SpectralFilterBankScaleType, WaveletContinueType
from audioflux.utils import note_to_hz

# Get a 220Hz's audio file path
sample_path = af.utils.sample_path('220')

# Read audio data and sample rate
audio_arr, sr = af.read(sample_path)
audio_arr = audio_arr[:4096]

cwt_obj = af.CWT(num=84, radix2_exp=12, samplate=sr, low_fre=note_to_hz('C1'),
                 bin_per_octave=12, wavelet_type=WaveletContinueType.MORSE,
                 scale_type=SpectralFilterBankScaleType.OCTAVE)

cwt_spec_arr = cwt_obj.cwt(audio_arr)

synsq_obj = af.Synsq(num=cwt_obj.num,
                     radix2_exp=cwt_obj.radix2_exp,
                     samplate=cwt_obj.samplate)

synsq_arr = synsq_obj.synsq(cwt_spec_arr,
                            filter_bank_type=cwt_obj.scale_type,
                            fre_arr=cwt_obj.get_fre_band_arr())

# Show CWT
fig, ax = plt.subplots(figsize=(7,4))
img = fill_spec(np.abs(cwt_spec_arr), axes=ax,
                x_coords=cwt_obj.x_coords(),
                y_coords=cwt_obj.y_coords(),
                x_axis='time', y_axis='log',
                title='CWT')
fig.colorbar(img, ax=ax)
# Show Synsq
fig, ax = plt.subplots(figsize=(7,4))
img = fill_spec(np.abs(synsq_arr), axes=ax,
                x_coords=cwt_obj.x_coords(),
                y_coords=cwt_obj.y_coords(),
                x_axis='time', y_axis='log',
                title='Synsq')
fig.colorbar(img, ax=ax)

plt.show()
```

<img src='image/demo_cwt.png'  width="415"  /><img src='image/demo_synsq.png'  width="415"  />


### Other examples

- [CQT & Chroma](docs/examples.md#cqt--chroma)
- [Different Wavelet Type](docs/examples.md#different-wavelet-type)
- [Spectral Features](docs/examples.md#spectral-features)
- [Pitch Estimate](docs/examples.md#pitch-estimate)
- [Onset Detection](docs/examples.md#onset-detection)
- [Harmonic Percussive Source Separation](docs/examples.md#harmonic-percussive-source-separation)

More example scripts are provided in the [Documentation](https://audioflux.top/) section.

## Installation
![language](https://img.shields.io/badge/platform-iOS%20|%20android%20|%20macOS%20|%20linux%20|%20windows%20-lyellow.svg) 

The library is cross-platform and currently supports Linux, macOS, Windows, iOS and Android systems.

### Python Package Intsall  

Using PyPI: 

```
$ pip install audioflux 
```

<!--Using Anaconda: 

```
$ conda install -c conda-forge audioflux
```-->


<!--Read installation instructions:
https://audioflux.top/install-->


### iOS build

To compile iOS on a Mac, Xcode Command Line Tools must exist in the system:  

- Install the full Xcode package
- install Xcode Command Line Tools when triggered by a command or run xcode-select command: 

```
$ xcode-select --install 
```
Enter the **`audioFlux`** project **`scripts`** directory and switch to the current directory, run the following script to build and compile:

```
$ ./build_iOS.sh
```

Build  and compile successfully, the project build compilation results are in the **`build`** folder

### Android build

The current system development environment needs to be installed [**android NDK**](https://developer.android.com/ndk), ndk version>=16,after installation, set the environment variable ndk path.  
 
For example, ndk installation path is `~/Android/android-ndk-r16b`:  

```
$ export NDK_ROOT=~/Android/android-ndk-r16b
$ export PATH=$NDK_ROOT:$PATH
```

Android **`audioFlux`** build uses [**fftw**](https://www.fftw.org/) library to accelerate performance, compile the single-floating point version for android platform. fftw lib successful compilation, copy to  **`audioFlux`** project **`scripts/android/fftw3`** directory.

Enter the **`audioFlux`** project **`scripts`** directory and switch to the current directory, run the following script to build and compile:

```
$ ./build_android.sh
```

Build  and compile successfully, the project build compilation results are in the **`build`** folder


### Building from source

For Linux, macOS, Windows systems. Read installation instructions:

* [docs/installing.md](docs/installing.md)

## Documentation

Documentation of the package can be found online:

[https://audioflux.top](https://audioflux.top/)

## Contributing
We are more than happy to collaborate and receive your contributions to **`audioFlux`**. If you want to contribute, please fork the latest git repository and create a feature branch. Submitted requests should pass all continuous integration tests.

You are also more than welcome to suggest any improvements, including proposals for need help, find a bug, have a feature request, ask a general question, new algorithms. <a href="https://github.com/libAudioFlux/audioFlux/issues/new">Open an issue</a>

## Citing

If you want to cite **`audioFlux`** in a scholarly work, there are two ways to do it.

- If you are using the library for your work, for the sake of reproducibility, please cite
  the version you used as indexed at Zenodo:

    [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7548288.svg)](https://doi.org/10.5281/zenodo.7548288)

## License
audioFlux project is available MIT License.


