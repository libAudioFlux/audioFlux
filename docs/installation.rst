Installation
============

The library is cross-platform and currently supports Linux, macOS, Windows, iOS and Android systems.

Python Package Install
----------------------

To install the ``audioFlux`` package, you must have Python 3.6 or newer

We recommend installing using the released python package:

Using PyPI
^^^^^^^^^^
 ::

    $ pip install audioflux


Using conda
^^^^^^^^^^^
 ::

    $ conda install -c tanky25 -c conda-forge audioflux


Building for mobile
-------------------

iOS build
^^^^^^^^^

To compile iOS on a Mac, Xcode Command Line Tools must exist in the system:

- Install the full Xcode package
- Install Xcode Command Line Tools when triggered by a command or run xcode-select command::

    $ xcode-select --install


Enter the ``audioFlux`` project ``scripts`` directory and switch to the current directory, run the following script to build and compile::

    $ ./build_iOS.sh


Build and compile successfully, the project build compilation results are in the ``build`` folder

Android build
^^^^^^^^^^^^^
The current system development environment needs to be installed `android NDK <https://developer.android.com/ndk>`_, ndk version>=16,after installation, set the environment variable ndk path.

For example, ndk installation path is ``~/Android/android-ndk-r16b``::

    $ export NDK_ROOT=~/Android/android-ndk-r16b
    $ export PATH=$NDK_ROOT:$PATH


Android ``audioFlux`` build uses `fftw <https://www.fftw.org/>`_ library to accelerate performance, compile the single-floating point version for android platform. fftw lib successful compilation, copy to  ``audioFlux`` project ``scripts/android/fftw3`` directory.

Enter the ``audioFlux`` project ``scripts`` directory and switch to the current directory, run the following script to build and compile::

    $ ./build_android.sh


Build and compile successfully, the project build compilation results are in the ``build`` folder


Building from source
--------------------

1. **Prerequisites**

    Install fft library

    * Linux:

        Install **fftw-3f**, detail see: `Install fftw <https://www.fftw.org/download.html>`_

            For Ubuntu::

                sudo apt install libfftw3-dev


            For CentOS::

                sudo yum install fftw3-devel


    * macOS:

        If version>=13, need to install Xcode and ``xcode-select``.

2. **Build and install**::

    $ python setup.py build
    $ python setup.py install


