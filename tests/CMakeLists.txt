# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
message("\n Adding the test scripts")

# Add the files to be compiled into a library
aux_source_directory(.					SOURCE_FILES)
message(STATUS						"${SOURCE_FILES}")

# Add executables
add_executable(			test_quality		TestEstimatorQuality.cc)
target_link_libraries(		test_quality		Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			test_convergence	TestEstimatorConvergence.cc)
target_link_libraries(		test_convergence	Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			test_alphacomplex	TestAlphaComplex.cc)
target_link_libraries(		test_alphacomplex	Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			test_optimal		TestOptimalParameters.cc)
target_link_libraries(		test_optimal		Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			test_hull		TestConvexeHull.cc)
target_link_libraries(		test_hull		Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			test_amplitude		TestAmplitudeRecon.cc)
target_link_libraries(		test_amplitude		Synapse
							${ROOT_LIBRARIES}
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			sample_string		SampleFromString.cc)
target_link_libraries(		sample_string		${ROOT_LIBRARIES}
							-lGeom)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -O3")
