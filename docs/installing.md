# Installation

![language](https://img.shields.io/badge/platform-%20Linux%20%7C%20macOS%20%7C%20Windows%20%7C%20iOS%20%7C%20Android%20-lyellow.svg)

The library is cross-platform and currently supports Linux, macOS, Windows, iOS and Android systems.

- [Python Package Install](#python-package-install)
- [iOS build](#ios-build)
- [Android build](#android-build)
- [Building from source](#building-from-source)

## Python Package Install

To install the **audioFlux** package, Python >=3.6, using the released python package.

Using PyPI:

```
$ pip install audioflux 
```

Using Anaconda:

```
$ conda install -c tanky25 -c conda-forge audioflux
```

## iOS build

To compile iOS on a Mac, Xcode Command Line Tools must exist in the system:

- Install the full Xcode package
- install Xcode Command Line Tools when triggered by a command or run xcode-select command:

```
$ xcode-select --install 
```

Enter the **`audioFlux`** project **`scripts`** directory and switch to the current directory, run the following script
to build and compile:

```
$ ./build_iOS.sh
```

Build and compile successfully, the project build compilation results are in the **`build`** folder

## Android build

The current system development environment needs to be installed [**android NDK**](https://developer.android.com/ndk),
ndk version>=16,after installation, set the environment variable ndk path.

For example, ndk installation path is `~/Android/android-ndk-r16b`:

```
$ export NDK_ROOT=~/Android/android-ndk-r16b
$ export PATH=$NDK_ROOT:$PATH
```

Android **`audioFlux`** build uses [**fftw**](https://www.fftw.org/) library to accelerate performance, compile the
single-floating point version for android platform. fftw lib successful compilation, copy to  **`audioFlux`**
project **`scripts/android/fftw3`** directory.

Enter the **`audioFlux`** project **`scripts`** directory and switch to the current directory, run the following script
to build and compile:

```
$ ./build_android.sh
```

Build and compile successfully, the project build compilation results are in the **`build`** folder

## Building from source

#### Linux build

1. Install compilation dependencies

    ```shell
    # For Ubuntu:
    $ sudo apt install cmake clang
    
    # For CentOS:
    $ sudo yum install -y cmake clang
    ```

2. Installing **MKL** lib dependencies on Linux

   You can
   use [this installation document](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux)
   to install MKL

   After installation, you need to set the environment variable `MKL_INCLUDE_PATH` and `MKL_LIB_PATH` for MKL.

   For example:

   ```shell
   $ export MKL_INCLUDE_PATH=/opt/intel/oneapi/mkl/latest/include
   $ export MKL_LIB_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
   ```
<!-- 

* Install MKL using pip.
   ```shell
   $ pip3 install mkl mkl-include
   ```

* (Optional) If you use pip in a non-system directory, such as a virtual environment created by virtualenv, you need to
  set the environment variables `MKL_INCLUDE_PATH` and `MKL_LIB_PATH` of the MKL directory
   ```shell
   $ export MKL_INCLUDE_PATH=/opt/intel/oneapi/mkl/latest/include
   $ export MKL_LIB_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
   ```

* Add soft link for compilation.
-->

3. Python setup:

    ```shell
    $ python setup.py build
    $ python setup.py install
    ```

#### macOS build

1. Installing dependencies on macOS

    * Install Command Line Tools for Xcode. Even if you install Xcode from the app store you must configure command-line
      compilation by running:

        ```shell
        $ xcode-select --install
        ```

    * Install **llvm**.
        ```shell
        $ brew install llvm
      
        $ # After installation, you need to set the environment variable C_COMPILER_PATH and CXX_COMPILER_PATH for MKL. For example:
        $ export C_COMPILER_PATH=/usr/local/opt/llvm/bin/clang
        $ export CXX_COMPILER_PATH=/usr/local/opt/llvm/bin/clang++
        ```

2. Python setup:

   ```shell
   $ python setup.py build
   $ python setup.py install
   ```

#### Windows build

Building from source is currently not supported. Only supports pip installation. If you need to build from source, The
easiest way to build audioFlux is by cross-compilation on Linux/macOS using MinGW.

 ```shell
 $ python setup.py build_py_win
 $ python setup.py install
 ```
