#! /bin/sh

cd ../
rm -r build
mkdir build && cd build

cp ../scripts/windows/fftw3/fftw3.h ../src/
cp ../scripts/windows/fftw3/libfftw3f-3.dll ../src/

cmake ../src -BwindowBuild \
	-DCMAKE_TOOLCHAIN_FILE=../scripts/win64.cmake \
	-DCMAKE_SYSTEM_NAME=win64 \

cmake --build windowBuild --config Release 

rm ../src/fftw3.h
rm ../src/libfftw3f-3.dll