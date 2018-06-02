# Check if cmake has the required version
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
message("\n Adding the main executables")

# Add the files to be compiled
aux_source_directory(.					SOURCE_FILES)
message(STATUS						"${SOURCE_FILES}")

# Add executables
if (MAUS_FOUND)
  add_executable(		import_data		ImportMAUSData.cc)
  target_link_libraries(	import_data		Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
        	                              		${MAUS_LIBRARIES})

  add_executable(		import_sim		ImportMAUSSimulation.cc)
  target_link_libraries(	import_sim		Synapse
							${ROOT_LIBRARIES}
                                      			${MAUS_LIBRARIES})

  add_executable(		import_minimal		ImportMAUSMinimal.cc)
  target_link_libraries(	import_minimal		Synapse
							${ROOT_LIBRARIES}
                                      			${MAUS_LIBRARIES})

  add_executable(		convert_to_json		ConvertToJson.cc)
  target_link_libraries(	convert_to_json		Synapse
							${ROOT_LIBRARIES})
endif (MAUS_FOUND)

add_executable(			import_g4bl		ImportG4BL.cc)
target_link_libraries(		import_g4bl		Synapse
							${ROOT_LIBRARIES})

add_executable(			phase_space		PhaseSpace.cc)
target_link_libraries(		phase_space		Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			reweight		ReweightBeam.cc)
target_link_libraries(		reweight		Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			syst			SystematicCorrections.cc)
target_link_libraries(		syst			Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			amplitudes		Amplitudes.cc)
target_link_libraries(		amplitudes		Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			density			DensityProfile.cc)
target_link_libraries(		density			Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			profiles		Profiles.cc)
target_link_libraries(		profiles		Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			toy_sim			ToySimulation.cc)
target_link_libraries(		toy_sim			Synapse
							${ROOT_LIBRARIES}
							-lSpectrum
							-lMinuit
							-lGeom
							qhull_r
							qhullcpp)

add_executable(			compare_sets		CompareSets.cc)
target_link_libraries(		compare_sets		Synapse
							${ROOT_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -O3")
