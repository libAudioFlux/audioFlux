cmake_minimum_required(VERSION 3.5)
project(audioflux)

aux_source_directory(vector vector)
aux_source_directory(util util)
aux_source_directory(dsp dsp)
aux_source_directory(classic classic)
aux_source_directory(filterbank filterbank)
aux_source_directory(feature feature)
aux_source_directory(mir mir)
aux_source_directory(. SOURCES)

list(APPEND SOURCES ${vector} ${util} ${dsp} ${classic} ${filterbank} ${feature} ${mir})

if(${CMAKE_SYSTEM_NAME} MATCHES "iOS")
	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_ACCELERATE")

	# set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")

	add_library(audioflux
				${SOURCES})

	target_link_libraries(audioflux
						"-framework Accelerate")

	install(TARGETS audioflux DESTINATION lib)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "darwin")
	set(CMAKE_C_COMPILER $ENV{C_COMPILER_PATH})
	set(CMAKE_CXX_COMPILER $ENV{CXX_COMPILER_PATH})
	if (DEFINED ENV{AF_BUILD_PY_BDIST})
		message("Build py_bdist version")

		set(ompPath ../scripts/macOS/openmp)
		include_directories(${ompPath})
		link_directories(${ompPath})
	endif()

	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -O3 -DHAVE_ACCELERATE -DHAVE_OMP -fopenmp ")

	add_library(audioflux
				SHARED
				${SOURCES})

	target_link_libraries(audioflux
						"-framework Accelerate"
						omp
						)

	install(TARGETS audioflux DESTINATION lib)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "android")
	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c99  -DHAVE_FFTW3F ")

	# strip
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")

	if(${CMAKE_ANDROID_ARCH_ABI} MATCHES "armeabi-v7a")
		set(fftwPath ../scripts/android/fftw3/armeabi-v7a)
	elseif(${CMAKE_ANDROID_ARCH_ABI} MATCHES "arm64-v8a")
		set(fftwPath ../scripts/android/fftw3/arm64-v8a)
	elseif(${CMAKE_ANDROID_ARCH_ABI} MATCHES "armeabi")
		set(fftwPath ../scripts/android/fftw3/armeabi)
	endif()

	include_directories(${fftwPath})
	link_directories(${fftwPath})

	# add dylib, include source path
	add_library(audioflux
				SHARED
				${SOURCES})

	# search system lib, ndk-10e is not support
	find_library( # Sets the name of the path variable.
				syslibs
				log
				android)

	# link target lib
	target_link_libraries(audioflux
						 fftw3f
						 log)

	install(TARGETS audioflux DESTINATION lib)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "linux")
	set(CMAKE_C_COMPILER "/usr/bin/clang")
	set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
	if (DEFINED ENV{AF_BUILD_PY_BDIST})
		message("Build py_bdist version")
		set(OMP_PATH "../scripts/linux/openmp")
		set(MKL_INCLUDE_PATH "../scripts/linux/mkl")
		set(MKL_LIB_PATH "../scripts/linux/mkl")

		include_directories(${OMP_PATH})
		link_directories(${OMP_PATH})
	else()
		set(MKL_INCLUDE_PATH $ENV{MKL_INCLUDE_PATH})
		set(MKL_LIB_PATH $ENV{MKL_LIB_PATH})
		message("Using MKL_INCLUDE_PATH: ${MKL_INCLUDE_PATH}")
		message("Using MKL_LIB_PATH: ${MKL_LIB_PATH}")
	endif()

	if (MKL_LIB_PATH)
		set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -O3 -DHAVE_MKL -DHAVE_OMP -fopenmp=libomp -Wl,-rpath,${MKL_LIB_PATH} ")
	else()
		set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -O3 -DHAVE_MKL -DHAVE_OMP -fopenmp=libomp ")
	endif()

	if (MKL_INCLUDE_PATH)
		include_directories(${MKL_INCLUDE_PATH})
	endif()
	if (MKL_LIB_PATH)
		link_directories(${MKL_LIB_PATH})
	endif()
	include_directories(${PROJECT_SOURCE_DIR})
	link_directories(${PROJECT_SOURCE_DIR})

	# add so, include source path
	add_library(audioflux
				SHARED
				${SOURCES})

	# link target lib
	target_link_libraries(audioflux
						 mkl_intel_lp64
						 mkl_core
						 mkl_intel_thread
						 omp
						 m)

	install(TARGETS audioflux DESTINATION lib)

elseif (${CMAKE_SYSTEM_NAME} MATCHES "win64")
	# -Wl,-rpath=.
	# -Wl,--out-implib,./audioflux.lib  MSVC lib
	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_FFTW3F -Wl,-rpath=.  ")

	set(BUILD_SHARED_LIBS true)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

	set(fftwPath ../scripts/windows/fftw3)
	include_directories(${fftwPath})

	link_directories(${PROJECT_SOURCE_DIR})
	include_directories(${PROJECT_SOURCE_DIR})

	# add so, include source path
	add_library(audioflux
				SHARED
				${SOURCES})

	# link target lib
	target_link_libraries(audioflux
						 fftw3f-3
						 )

	install(TARGETS audioflux DESTINATION lib)

endif()






