

set(CMAKE_SYSTEM_NAME win64)

set(tools /usr/local/Cellar/mingw-w64/10.0.0_3/toolchain-x86_64)
set(CMAKE_C_COMPILER ${tools}/bin/x86_64-w64-mingw32-gcc)
set(CMAKE_CXX_COMPILER ${tools}/bin/x86_64-w64-mingw32-g++)


set(BUILD_SHARED_LIBS true)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)