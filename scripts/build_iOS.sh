#! /bin/sh

cd ../
rm -r build
mkdir build && cd build 

# armv7;armv7s;arm64;i386;x86_64

cmake ../src -BiOSBuild -GXcode \
    -DCMAKE_SYSTEM_NAME=iOS \
    "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64" \
    -DCMAKE_OSX_DEPLOYMENT_TARGET=11.0 \
    -DCMAKE_INSTALL_PREFIX=`pwd`/_install \
    -DCMAKE_XCODE_ATTRIBUTE_ONLY_ACTIVE_ARCH=NO \
    -DCMAKE_IOS_INSTALL_COMBINED=YES \

cmake --build iOSBuild --config Release 
# cmake --build iOSBuild --config Debug 