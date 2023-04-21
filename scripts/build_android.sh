#! /bin/sh

cd ../
rm -r build
mkdir build && cd build 

# -DANDROID new way, must have toolchain
# -DCMAKE_ cmake way, compatibility better

# -DCMAKE_TOOLCHAIN_FILE=$NDK_ROOT/build/cmake/android.toolchain.cmake \
# -DBUILD_SHARED_LIBS=1 \
# -DCMAKE_BUILD_TYPE=Release \
# -DANDROID_ABI=arm64-v8a armeabi armeabi-v7a \

# -DANDROID_ARM_NEON=TRUE	 \
# -DCMAKE_ANDROID_STL_TYPE=gnustl_static \
# -DCMAKE_ANDROID_STL_TYPE=c++_shared \

# -DANDROID_NATIVE_API_LEVEL=21 \
# -DCMAKE_ANDROID_API=21
# -DANDROID_TOOLCHAIN=clang \
# -DANDROID_NDK=$NDK_ROOT \

cmake ../src -Barm64-v8a \
	-DCMAKE_ANDROID_API=21 \
	-DCMAKE_SYSTEM_NAME=android \
	"-DCMAKE_ANDROID_ARCH_ABI=arm64-v8a"  \
	-DCMAKE_ANDROID_NDK=$NDK_ROOT \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_ANDROID_NDK_TOOLCHAIN_VERSION=clang \
	-DCMAKE_ANDROID_SKIP_ANT_STEP=1 \

# make -j8
cmake --build arm64-v8a --config Release 

cmake ../src -Barmeabi-v7a \
	-DCMAKE_ANDROID_API=21 \
	-DCMAKE_SYSTEM_NAME=android \
	"-DCMAKE_ANDROID_ARCH_ABI=armeabi-v7a"  \
	-DCMAKE_ANDROID_NDK=$NDK_ROOT \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_ANDROID_NDK_TOOLCHAIN_VERSION=clang \
	-DCMAKE_ANDROID_SKIP_ANT_STEP=1 \


# make -j8
cmake --build armeabi-v7a --config Release

cmake ../src -Barmeabi \
	-DCMAKE_ANDROID_API=21 \
	-DCMAKE_SYSTEM_NAME=android \
	"-DCMAKE_ANDROID_ARCH_ABI=armeabi"  \
	-DCMAKE_ANDROID_NDK=$NDK_ROOT \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_ANDROID_NDK_TOOLCHAIN_VERSION=clang \
	-DCMAKE_ANDROID_SKIP_ANT_STEP=1 \


# make -j8
cmake --build armeabi --config Release











	

