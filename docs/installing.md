
## Compiling from source

#### Linux build

1. Installing dependencies on Linux

	For Ubuntu: 
		
	```
	sudo apt install libfftw3-dev
	```
		
	For CentOS: 
		
	```
	sudo yum install -y fftw3-devel
	```
	
2. Python setup:

	```
	$ python setup.py build
	$ python setup.py install
	```
	
#### macOS build

1. Installing dependencies on macOS

	Install Command Line Tools for Xcode. Even if you install Xcode from the app store you must configure command-line compilation by running:
	
	```
	xcode-select --install
	```

2. Python setup:

	```
	$ python setup.py build
	$ python setup.py install
	```
 
#### Windows build

Building from source is currently not supported. Only supports pip installation. If you need to build from source, The easiest way to build audioFlux is by cross-compilation on Linux using MinGW.
