cmake_minimum_required(VERSION 3.5)
project(FluxLib)

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
	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_ACCELERATE")

	add_library(audioflux
				SHARED
				${SOURCES})

	target_link_libraries(audioflux
						"-framework Accelerate")

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
	set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DHAVE_FFTW3F")

	# add so, include source path
	add_library(audioflux
				SHARED
				${SOURCES})

	# link target lib
	target_link_libraries(audioflux
						 fftw3f
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






