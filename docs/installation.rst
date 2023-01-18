Installation
============

The library is cross-platform and currently supports Linux, macOS, Windows, iOS and Android systems.

Python Package Install
----------------------

To install the ``audioFlux`` package, you must have Python 3.6 or newer

Using PyPI::

    $ pip install audioflux

Using Anaconda::

    $ conda install -c conda-forge audioflux


Building from source::

    $ python setup.py build
    $ python setup.py install

.. note::
    Installation of fft library

    * Linux:

        Install **fftw-3.x**

        For Ubuntu::

            sudo apt install libfftw3-dev

        For CentOS::

            sudo yum install -y fftw3-devel

        For FreeBSD::

            sudo pkg install fftw3

        Detail See: `Install fftw <https://www.fftw.org/download.html>`_

    * Windows:

        Install **fftw-3.x**

        Detail See: `Install Windows fftw <https://www.fftw.org/install/windows.html>`_

    * macOS:

        If version>=13, need to install Xcode and ``xcode-select``.

iOS build
---------

To compile iOS on a Mac, Xcode Command Line Tools must exist in the system:

- Install the full Xcode package
- install Xcode Command Line Tools when triggered by a command or run xcode-select command::

    $ xcode-select --install

Enter the ``audioFlux`` project ``scripts`` directory and switch to the current directory,
run the following script to build and compile::

    $ ./buildFluxLib_iOS.sh

If the system has permission problems, run::

    $ sudo ./buildFluxLib_iOS.sh


Build  and compile successfully, the project build compilation results are in the ``build`` folder

Android build
-------------

The current system development environment needs to be installed `android NDK <https://developer.android.com/ndk>`_,
ndk version>=16,after installation, set the environment variable ndk path.

For example, ndk installation path is `~/Android/android-ndk-r16b`::

    $ NDK_ROOT=~/Android/android-ndk-r16b
    $ export PATH=$NDK_ROOT:$PATH


Android ``audioFlux`` build uses `fftw <https://www.fftw.org/>`_ library to accelerate performance, compile
the single-floating point version for android platform. fftw lib successful compilation, copy to
``audioFlux`` project ``scripts/android/fftw3`` directory.

Enter the ``audioFlux`` project ``scripts`` directory and switch to the current directory, run the following
script to build and compile::

    $ ./buildFluxLib_android.sh


Build  and compile successfully, the project build compilation results are in the ``build`` folder
