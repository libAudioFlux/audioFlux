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

Linux build
^^^^^^^^^^^

1. Install compilation dependencies::

    # For Ubuntu:
    $ sudo apt install cmake clang

    # For CentOS:
    $ sudo yum install -y cmake clang

2. Installing **MKL** lib dependencies on Linux (Optional, Recommended installation)

   You can
   use `this installation document <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux>`_
   to install MKL

   After installation, you need to set the environment variable `MKL_INCLUDE_PATH` and `MKL_LIB_PATH` for MKL.

   For example::

   $ export MKL_INCLUDE_PATH=/opt/intel/oneapi/mkl/latest/include
   $ export MKL_LIB_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64



3. Python setup::

    $ python setup.py build
    $ python setup.py install

macOS build
^^^^^^^^^^^

1. Installing dependencies on macOS

    * Install Command Line Tools for Xcode. Even if you install Xcode from the app store you must configure command-line
      compilation by running::

        $ xcode-select --install


    * Install **llvm**::

        $ brew install llvm

        $ # After installation, you need to set the environment variable C_COMPILER_PATH and CXX_COMPILER_PATH for MKL. For example:
        $ export C_COMPILER_PATH=/usr/local/opt/llvm/bin/clang
        $ export CXX_COMPILER_PATH=/usr/local/opt/llvm/bin/clang++


2. Python setup::

   $ python setup.py build
   $ python setup.py install


Windows build
^^^^^^^^^^^^^

Building from source code is currently not supported on Windows; it is recommended to use pip installation.

If you need to build from source, The easiest way to build audioFlux is by cross-compilation on Linux/macOS using MinGW.::

 $ python setup.py build_py_win
 $ python setup.py install

