#! /bin/sh

cd ../
rm -r build
mkdir build && cd build 

# -GXcode
# "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64" \
#  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13 \
# -DCMAKE_INSTALL_PREFIX=`pwd`/_install \
# -DCMAKE_XCODE_ATTRIBUTE_ONLY_ACTIVE_ARCH=NO \

cmake ../src -BmacOSBuild  \
	-DCMAKE_SYSTEM_NAME=darwin \
#	"-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64" \


cmake --build macOSBuild --config Release