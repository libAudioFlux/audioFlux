
## Building for mobile

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
