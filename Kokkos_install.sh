export PREFIX=${CONDA_PREFIX}
export WORKDIR=${PWD}
export CABANASRC=${CONDA_PREFIX}/cabana
export KOKKOSSRC=${CONDA_PREFIX}/kokkos

conda env config vars set PREFIX=${PREFIX} WORKDIR=${WORKDIR} CPATH=${PREFIX}/include CABANASRC=${CABANASRC} KOKKOSSRC=${KOKKOSSRC}

git clone https://github.com/kokkos/kokkos.git $KOKKOSSRC
cd $KOKKOSSRC
git checkout 4.1.00
mkdir -p $KOKKOSSRC/build && cd $KOKKOSSRC/build

cmake -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_CXX_COMPILER=g++ \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_EXTENSIONS=Off \
      -DKokkos_ENABLE_COMPILER_WARNINGS=ON \
      -DKokkos_ENABLE_CUDA=Off \
      -DKokkos_ENABLE_CUDA_LAMBDA=Off \
      -DKokkos_ENABLE_OPENMP=On \
      -DKokkos_ENABLE_SERIAL=On \
      -DKokkos_ENABLE_TESTS=Off \
      -DKokkos_ARCH_AMPERE80=Off $KOKKOSSRC
cmake --build . --target install -j

git clone https://github.com/ECP-copa/Cabana.git $CABANASRC
cd $CABANASRC  
git checkout 2b642f8dbba760e2008f8b81a4196e6adf100f42
mkdir -p $CABANASRC/build && cd $CABANASRC/build
    cmake -D CMAKE_BUILD_TYPE="Release" \
          -D CMAKE_CXX_STANDARD=17 \
          -D CMAKE_PREFIX_PATH=$CABANASRC \
          -D CMAKE_INSTALL_PREFIX=$PREFIX \
          -D CMAKE_CXX_COMPILER=g++\
          -D Cabana_REQUIRE_CUDA=OFF \
          -D Cabana_ENABLE_TESTING=OFF \
          -D Cabana_ENABLE_EXAMPLES=OFF $CABANASRC
cmake --build . --target install -j
