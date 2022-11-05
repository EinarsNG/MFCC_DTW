project(MFFC_DTW)

set(CMAKE_BINARY_DIR bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_compile_options(-g)

set (UTILS
	${CMAKE_SOURCE_DIR}/src/utils.cpp
	${CMAKE_SOURCE_DIR}/src/steps.cpp
)

add_executable(preemph
	${CMAKE_SOURCE_DIR}/steps/preemphasis.cpp
	${UTILS}
)

add_executable(framing
	${CMAKE_SOURCE_DIR}/steps/framing.cpp
	${UTILS}
)

add_executable(windowing
	${CMAKE_SOURCE_DIR}/steps/windowing.cpp
	${UTILS}
)
add_executable(dft
	${CMAKE_SOURCE_DIR}/steps/dft.cpp
	${UTILS}
)

add_executable(powspec
	${CMAKE_SOURCE_DIR}/steps/powerspectrum.cpp
	${UTILS}
)

add_executable(filterbanks
	${CMAKE_SOURCE_DIR}/steps/filterbanks.cpp
	${UTILS}
)

add_executable(log
	${CMAKE_SOURCE_DIR}/steps/logarithm.cpp
	${UTILS}
)

add_executable(dct
	${CMAKE_SOURCE_DIR}/steps/dct.cpp
	${UTILS}
)

add_executable(dtw
	${CMAKE_SOURCE_DIR}/steps/dtw.cpp
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
target_link_libraries(dtw fftw3)
