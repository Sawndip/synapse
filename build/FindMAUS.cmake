# - Try to find MAUS
# Once done this will define
#
#  MAUS_FOUND - system has MAUS
#  MONITOR_FOUND - system has MDMonitor
#  MAUS_INCLUDE_DIR - the MAUS include directory
#  MAUS_LIBRARIES - The libraries needed to use MAUS
#  MAUS_DEFINITIONS - Compiler switches required for using MAUS
#

MESSAGE("\n Looking for MAUS...")
                   # Find MAUS_DIR
#message("MAUS_LIBRARIES : ${MAUS_LIBRARIES}")


if(NOT EXISTS $ENV{MAUS_ROOT_DIR})

  message(STATUS "MAUS NOT FOUND!  $MAUS_ROOT_DIR environmental variable is not set properly.
Did you try running: source env.sh? \nThe MAUS importer will not be built.")
  return()

endif(NOT EXISTS $ENV{MAUS_ROOT_DIR})


FIND_PATH(MAUS_LIBRARY_DIR NAMES    libMausCpp.so          PATHS
                                                           $ENV{MAUS_ROOT_DIR}/build/
                                                           NO_DEFAULT_PATH)

if (MAUS_LIBRARY_DIR)

  set(MAUS_FOUND TRUE)
  set(MAUS_LIBRARY_DIR  ${MAUS_LIBRARY_DIR}
                        $ENV{MAUS_THIRD_PARTY}/third_party/install/lib
                        $ENV{MAUS_THIRD_PARTY}/third_party/install/lib64
                       )


  set(MAUS_INCLUDE_DIR  $ENV{MAUS_THIRD_PARTY}/third_party/install/include
                        $ENV{MAUS_THIRD_PARTY}/third_party/install/include/python2.7
			$ENV{MAUS_ROOT_DIR}
                        $ENV{MAUS_ROOT_DIR}/src/common_cpp
                        $ENV{MAUS_ROOT_DIR}/src/legacy
                        $ENV{MAUS_ROOT_DIR}/src/map
                       )

  set(MAUS_LIBRARIES    MausCpp
                        InputCppDAQOfflineData
                        MapCppTOFDigits
                        MapCppTOFSlabHits
                        MapCppTOFSpacePoints
                        MapCppCkovDigits
                        MapCppKLDigits
                        MapCppKLCellHits
                        MapCppTrackerDigits
			MapCppTrackerClusterRecon
			MapCppTrackerSpacePointRecon
			MapCppTrackerPatternRecognition
			MapCppTrackerPRSeed
			MapCppTrackerTrackFit
			MapCppEMRPlaneHits
			MapCppEMRSpacePoints
                        MapCppEMRRecon
                        ReduceCppTofCalib
                        ReduceCppTOFPlot
                        )

  message(STATUS "found in  $ENV{MAUS_ROOT_DIR}")

  FIND_PATH(MONITOR_LIBRARY_DIR NAMES    libMDMonitor.a    PATHS
                                                           $ENV{MAUS_ROOT_DIR}/build/
                                                           NO_DEFAULT_PATH)

  if    (MONITOR_LIBRARY_DIR)

    set(MONITOR_FOUND    TRUE)
    set(MAUS_LIBRARIES   ${MAUS_LIBRARIES} MDMonitor)
    message(STATUS "  found Online Monitor")

  endif (MONITOR_LIBRARY_DIR)

#   message(STATUS "MAUS_LIBRARY_DIR:  ${MAUS_LIBRARY_DIR}")
#   message(STATUS "MAUS_INCLUDE_DIR:  ${MAUS_INCLUDE_DIR}")
#   message(STATUS "MAUS_LIBRARIES:  ${MAUS_LIBRARIES}")

else (MAUS_LIBRARY_DIR)

  message(FATAL_ERROR "Could NOT find MAUS!\n")

endif (MAUS_LIBRARY_DIR)


