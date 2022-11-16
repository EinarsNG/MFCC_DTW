cmake_minimum_required(VERSION 3.15)

project(MFCC_DTW)

include_directories(
	.
	${FFTW_INCLUDES}
)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(full)
add_subdirectory(steps)
