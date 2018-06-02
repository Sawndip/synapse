# Check if cmake has the required version
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
message("\n Creating the SynAPSE library...")

# Add the files to be compiled into a library
message("\n Creating the SynAPSE library...")
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/exc 		EXCEPTION_SRC)
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/utils	UTILS_SRC)
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/math 	MATH_SRC)
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/de 		DE_SRC)
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/beam 	BEAM_SRC)
aux_source_directory (	${CMAKE_BINARY_DIR}/../src/toy 		TOY_SRC)

list(APPEND FULL_SRC	${EXCEPTION_SRC}
    			${UTILS_SRC}
			${MATH_SRC}
			${DE_SRC}
			${BEAM_SRC}
			${TOY_SRC})

if (MAUS_FOUND)

  # If MAUS is found, include the MAUS importer tools
  aux_source_directory(	${CMAKE_BINARY_DIR}/../src/import 	IMPORT_SRC)
  list(APPEND FULL_SRC	${IMPORT_SRC})
endif (MAUS_FOUND)

message(STATUS 		"${FULL_SRC}")
add_library(Synapse 	${FULL_SRC})
