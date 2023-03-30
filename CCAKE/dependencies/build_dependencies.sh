#!/bin/bash

BASEDIR=$(git rev-parse --show-toplevel)

cd ${BASEDIR}
git submodule update --init --recursive


BUILD_DIR=${BASEDIR}/CCAKE/build
KOKKOS_BUILD_DIR=${BUILD_DIR}/kokkos
CABANA_BUILD_DIR=${BUILD_DIR}/cabana

export INSTALL_DIR=${BASEDIR}/local

#creates the build directories
mkdir -p $KOKKOS_BUILD_DIR
mkdir -p $CABANA_BUILD_DIR
mkdir -p $INSTALL_DIR

#Compile kokkos
KOKKOS_SRC_DIR="${BASEDIR}/CCAKE/dependencies/kokkos"
CABANA_SRC_DIR="${BASEDIR}/CCAKE/dependencies/cabana"

cd $KOKKOS_BUILD_DIR
if [ -f CMakeCache.txt ]; then rm CMakeCache.txt; fi
cmake \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DCMAKE_CXX_COMPILER=g++\
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_EXTENSIONS=Off \
  -DKokkos_ENABLE_COMPILER_WARNINGS=ON \
  -DKokkos_ENABLE_CUDA=Off \
  -DKokkos_ENABLE_CUDA_LAMBDA=Off \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ENABLE_SERIAL=On \
  -DKokkos_ENABLE_TESTS=Off \
  -DKokkos_ARCH_${ARCH}=Off \
  ${KOKKOS_SRC_DIR}

cmake --build . --target install -j

cd $CABANA_BUILD_DIR
if [ -f CMakeCache.txt ]; then rm CMakeCache.txt; fi
cmake \
   -D CMAKE_BUILD_TYPE="Debug" \
   -D CMAKE_PREFIX_PATH=$INSTALL_DIR \
   -D CMAKE_INSTALL_PREFIX=$INSTALL_DIR \
   -D CMAKE_CXX_COMPILER=g++\
   -D Cabana_REQUIRE_CUDA=OFF \
   -D Cabana_ENABLE_TESTING=OFF \
   -D Cabana_ENABLE_EXAMPLES=OFF \
   ${CABANA_SRC_DIR}

cmake --build . --target install -j

echo "Dependencies installed in ${INSTALL_DIR}"