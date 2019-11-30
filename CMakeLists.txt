# Check if cmake has the required version
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# Set name of the project
project(SynAPSE)

# Path to the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/../bin")

# Load some basic macros which alow to find the required modules
include(		${CMAKE_BINARY_DIR}/FindMAUS.cmake)
include(		${CMAKE_BINARY_DIR}/FindROOT.cmake)
include(		${CMAKE_BINARY_DIR}/BuildDictionary.cmake)

# C++ 11 flag
set(CMAKE_CXX_FLAGS 	"${CMAKE_CXX_FLAGS} -std=c++0x -Wall -O3")

# LIBRARY OUTPUT PATH
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/../lib")
include_directories(	${CMAKE_SOURCE_DIR})

# Add MAUS
include_directories(  	${MAUS_INCLUDE_DIR})
link_directories(     	${MAUS_LIBRARY_DIR})

# Add ROOT
include_directories(	${ROOT_INCLUDE_DIR})
link_directories   (	${ROOT_LIBRARY_DIR})

# Link the QHull library and include the source files
include_directories(   	${CMAKE_BINARY_DIR}/../third_party/qhull/src)
include_directories(   	${CMAKE_BINARY_DIR}/../third_party/qhull/src/qhull_r)
include_directories(   	${CMAKE_BINARY_DIR}/../third_party/qhull/src/libqhullcpp)
link_directories(     	${CMAKE_BINARY_DIR}/../third_party/qhull/lib)

# Include the source of NanoFLANN (KD tree indexation)
include_directories(	${CMAKE_BINARY_DIR}/../third_party/nanoflann/include)
include_directories(	${CMAKE_BINARY_DIR}/../third_party/nanoflann/examples)

# Build the main SynAPSE library, include the headers
include_directories(   	${CMAKE_BINARY_DIR}/../src/exc)
include_directories(   	${CMAKE_BINARY_DIR}/../src/utils)
include_directories(   	${CMAKE_BINARY_DIR}/../src/math)
include_directories(   	${CMAKE_BINARY_DIR}/../src/de)
include_directories(   	${CMAKE_BINARY_DIR}/../src/beam)
include_directories(   	${CMAKE_BINARY_DIR}/../src/toy)
if (MAUS_FOUND)
  include_directories(	${CMAKE_BINARY_DIR}/../src/import)
endif (MAUS_FOUND)
add_subdirectory(	${CMAKE_BINARY_DIR}/../src)

# Set the build flags
set(STOP_SIGNAL "${CMAKE_BINARY_DIR}/stop")
message(STATUS ${STOP_SIGNAL})

message(STATUS "C flags      : " ${CMAKE_C_FLAGS})
message(STATUS "CXX flags    : " ${CMAKE_CXX_FLAGS})
message(STATUS "Linker flags : " ${CMAKE_EXE_LINKER_FLAGS}\n)

# Build the test scripts
add_subdirectory(	${CMAKE_BINARY_DIR}/../tests)

# Build the main executables
add_subdirectory(	${CMAKE_BINARY_DIR}/../scripts)
