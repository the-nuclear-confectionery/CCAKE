#!/usr/bin/env bash

# Build containers
cd build-scripts
./build-base-img.sh
./build.sh
./build-dev.sh
mv *.sif ../images
cd ..

# Compiles release of CCAKE
mkdir -p ../CCAKE/build
if [ -f ../CCAKE/build/clean-cmake.sh ]; then
   cd ../CCAKE/build
   ./clean-cmake.sh
fi
cd ../../apptainer
singularity exec --bind ../CCAKE:/CCAKE images/ccake.sif bash -c "cd /CCAKE/build ; cmake .. ; make -j"
