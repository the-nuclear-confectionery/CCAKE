# CCAKE

# 1. Installation

For basic usage, the apptainer image can be used. It contains all the dependencies needed to run the code on CPU.
To build the image, run the following command inside its folder (skip this step if someone gave you the image ready):
``` bash
./build-base-img.sh
./build.sh
```
To run the image, do:
``` bash
./run.sh
```

If you want to build locally with GPU support, follow the instructions below.

## 1.1 Requirements

- Nvidia GPU (tested with RTX 3050 Laptop)
    * In theory, AMD and Intel GPUs should work as well with minimal
    adaptations.
- Nvidia HPC SDK 22.11 (More recent versions should work, but are not tested)
- CMake 3.4 or higher

## 1.2 Installing dependencies.

### 1.2.1 Download and install Nvidia HPC SDK 22.11 (Skip for installation on Delta)
- Choose a place where you have written permissions. In a local machine, `/opt` may be a good choice
(It may be necessary to use `sudo chown -R $USER:$USER /opt` to give yourself permissions to write there)
- Download the Nvidia HPC SDK 22.11 from https://developer.nvidia.com/hpc-sdk
- Unpack it in the chosen location with `tar -xvf <path_to_sdk>.tar.gz`
- Install the SDK with `./<path_to_sdk>/install`
- Add the following lines to your "~/.bashrc" at the end of the install:
``` bash
MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/man; export MANPATH
PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/bin:$PATH; export PATH
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/comm_libs/mpi/bin:$PATH
export MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/comm_libs/mpi/man
export CUDA_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/cuda/11.8
```

### 1.2.2 Download and install Kokkos
- Assuming you installed the SDK in `/opt`, create a "local" folder inside it. Kokkos will be installed there.
- Create the following env variables:
``` bash
KOKKOSSRC=<path/to/kokkos_src>
PREFIX=/opt/local
```
> **Note**: If you exit the terminal, you will need to set these variables again.
- Clone the Kokkos repository: `git clone git@github.com:kokkos/kokkos.git`
- Checkout to version 4.1.0: `git checkout 4.1.00`
- Create a build folder inside the Kokkos folder: `mkdir -p $KOKKOSSRC/build && cd $KOKKOSSRC/build`
- Edit `$KOKKOSSRC/bin/nvcc_wrapper` according to your architecture. Most likely, you will need to replace `sm_70` by `sm_80` or `sm_86` (according to the architecture being targeted).
- Configure kokkos from within build dir with 
```bash
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_CXX_COMPILER=$KOKKOSSRC/bin/nvcc_wrapper \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_EXTENSIONS=Off \
      -DCMAKE_BUILD_TYPE="Release" \
      -DKokkos_ENABLE_COMPILER_WARNINGS=ON \
      -DKokkos_ENABLE_CUDA=On \
      -DKokkos_ENABLE_CUDA_LAMBDA=On \
      -DKokkos_ENABLE_OPENMP=On \
      -DKokkos_ENABLE_SERIAL=On \
      -DKokkos_ENABLE_TESTS=Off \
      -DKokkos_ARCH_AMPERE86=On $KOKKOSSRC
```
> **Note**: Depending on your GPU architecture, you may need to change the `Kokkos_ARCH_AMPERE86` flag.
8. Build and install: `cmake --build . --target install -j`
### 1.2.3 Download and install Cabana
- Assuming you installed the SDK in `/opt`, create a "local" folder inside it. Kokkos will be installed there.
- Create the following env variables:
``` bash
CABANASRC=<Path/to/cabana/src>
PREFIX=/opt/local
```
> **Note**: If you exit the terminal, you will need to set these variables again. Also, if you are not on the same terminal used to
install kokkos, you will need to set the `KOKKOSSRC` variable again.
- Clone the Cabana repository: `git clone https://github.com/ECP-copa/Cabana.git`
- Create a build folder inside the Cabana folder: `mkdir -p $CABANASRC/build && cd $CABANASRC/build`
- Configure Cabana from within build dir with 
```bash
cmake -D CMAKE_BUILD_TYPE="Debug" \
      -D CMAKE_CXX_STANDARD=17 \
      -D CMAKE_PREFIX_PATH=$CABANASRC \
      -D CMAKE_INSTALL_PREFIX=$PREFIX \
      -D CMAKE_CXX_COMPILER=$KOKKOSSRC/bin/nvcc_wrapper \
      -D Cabana_REQUIRE_CUDA=ON \
      -D Cabana_ENABLE_TESTING=OFF \
      -D Cabana_ENABLE_EXAMPLES=OFF $CABANASRC    
```
- Build and install: `cmake --build . --target install -j`

**Note:**
If using Cabana 0.5.0, there may be linkage issues. In that case, use 
'''bash
if(NOT TARGET Cabana::cabanacore)
  find_package(Cabana)
endif()
'''
### 1.3 Build CCAKE
**Note:**
If installing on Delta, will have to explicity install gsl and yaml first.
- For gsl, use the command ' module load gsl/2.7.1'.
- For yaml, clone using 'git clone https://github.com/jbeder/yaml-cpp.git'. Then,
'''bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME ..
make
make install 
'''

- Create a build folder:  `mkdir CCAKE/build && cd CCAKE/build`
- Configure CCAKE: `cmake -DCMAKE_PREFIX_PATH=/opt/local -DCMAKE_BUILD_TYPE="Debug" .. && make -j`
