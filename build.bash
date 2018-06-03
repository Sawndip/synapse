#!/bin/bash

FILE_STD=install.log

if [ -z "$MAUS_ROOT_DIR" ]; then
   echo "WARNING: \$MAUS_ROOT_DIR is not set"
   echo "WARNING: the MAUS data structure importer will not be built"
fi

if [ -f $FILE_STD ];
then
    rm $FILE_STD
fi

echo
echo "   Welcome to the all-in-one SynAPSE installer script. "
echo "   You can get the details of the progress "
echo "   (or ensure it is doing something) by running:"
echo
echo "       tail -f $FILE_STD"
echo

while [[ $# > 1 ]]; do
  key="$1"

  case $key in
  -j|--num-threads)
    if expr "$2" : '-\?[0-9]\+$' >/dev/null; then
      NUM_THREADS="$2"
    fi
    shift
    ;;
  -v|--verbosity)
    if [ $2 -eq 0 ] || [ $2 -eq 1 ]; then
      BUILD_VERBOSITY="$2"
    fi
  shift
  ;;
  --use-maus-gcc)
    if [ "$2" = true ] || [ "$2" = false ]; then
      USE_MAUS_GCC="$2"
     fi
     shift
     ;;
  esac
  shift
done

# Set default to use the system GCC, export the CC and CXX is the MAUS GCC is requested
if [ -z "$USE_MAUS_GCC" ]; then
    USE_MAUS_GCC=false
fi
if [ "$USE_MAUS_GCC" != true ] && [ "$USE_MAUS_GCC" != false ]; then
    USE_MAUS_GCC=false
fi
if [ "$USE_MAUS_GCC" = true ]; then
  echo
  echo "Using MAUS GCC"
  if [ ! -z ${MAUS_ROOT_DIR} ]; then
      echo ${MAUS_ROOT_DIR}/third_party/install/bin/gcc
      export CC=${MAUS_ROOT_DIR}/third_party/install/bin/gcc
      echo ${MAUS_ROOT_DIR}/third_party/install/bin/g++
      export CXX=${MAUS_ROOT_DIR}/third_party/install/bin/g++
  else
      echo "FATAL: \$MAUS_ROOT_DIR is not set"
      echo "FATAL: cannot use the MAUS GCC, abort."
  fi
  echo
else
  echo
  echo "Using System GCC"
  echo
fi

# Set the default number of threads to 1
if [ -z "$NUM_THREADS" ]; then
  NUM_THREADS=1
fi
echo "Number of threads to build with:" "${NUM_THREADS}"
echo

# Set the default verbosity to 0
if [ -z "$BUILD_VERBOSITY" ]; then
  BUILD_VERBOSITY=0
fi

echo "Installing third party libraries locally..."
echo

# Build the QHull 2015.2 libraries if they are not found exist
if [ ! -d third_party/qhull ]; then
  pushd third_party/
  tar -xvf qhull*tgz
  mv qhull-2015.2 qhull
  pushd qhull
  make -j ${NUM_THREADS}
  popd

  echo
  echo "Built Qhull 2015.2 at third_party/qhull"
  echo
else
  echo "Qhull 2015.2 has been found at third_party/qhull"
  echo
fi

# Unpack the NanoFLANN 1.2.3 headers if they are not found
if [ ! -d third_party/nanoflann ]; then
  pushd third_party/
  tar -xvf nanoflann*tar.gz
  mv nanoflann-1.2.3 nanoflann
  popd

  echo
  echo "Unpacked NanoFLANN 1.2.3 at third_party/nanoflann"
  echo
else
  echo "NanoFLANN 1.2.3 has been found at third_party/nanoflann"
  echo
fi

# Build the SynAPSE code
echo "Will now build the SynAPSE source code to bin/"
echo

pushd build/
find . -not -name FindROOT.cmake -not -name FindMAUS.cmake -delete
cmake ..
make -j ${NUM_THREADS}
popd
