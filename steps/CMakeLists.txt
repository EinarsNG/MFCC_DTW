project(MFFC_DTW)

set(CMAKE_BINARY_DIR bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set (UTILS
	${CMAKE_SOURCE_DIR}/src/utils.cpp
	${CMAKE_SOURCE_DIR}/src/steps.cpp
)

add_executable(preemph
	${CMAKE_SOURCE_DIR}/src/preemphasis.cpp
	${UTILS}
)

add_executable(framing
	${CMAKE_SOURCE_DIR}/src/framing.cpp
	${UTILS}
)

add_executable(windowing
	${CMAKE_SOURCE_DIR}/src/windowing.cpp
	${UTILS}
)
add_executable(dft
	${CMAKE_SOURCE_DIR}/src/dft.cpp
	${UTILS}
)

add_executable(powspec
	${CMAKE_SOURCE_DIR}/src/powerspectrum.cpp
	${UTILS}
)

add_executable(filterbanks
	${CMAKE_SOURCE_DIR}/src/filterbanks.cpp
	${UTILS}
)

add_executable(log
	${CMAKE_SOURCE_DIR}/src/logarithm.cpp
	${UTILS}
)

add_executable(dct
	${CMAKE_SOURCE_DIR}/src/dct.cpp
	${UTILS}
)

target_link_libraries(preemph fftw3)
target_link_libraries(framing fftw3)
target_link_libraries(windowing fftw3)
target_link_libraries(dft fftw3)
target_link_libraries(powspec fftw3)
target_link_libraries(filterbanks fftw3)
target_link_libraries(log fftw3)
target_link_libraries(dct fftw3)
