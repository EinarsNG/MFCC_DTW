project(MFFC_DTW)

set(CMAKE_BINARY_DIR bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

set (UTILS
	${CMAKE_SOURCE_DIR}/src/utils.cpp
	${CMAKE_SOURCE_DIR}/src/steps.cpp
)

add_executable(compare
	${CMAKE_SOURCE_DIR}/full/compare.cpp
	${UTILS}
)

target_link_libraries(compare fftw3)
