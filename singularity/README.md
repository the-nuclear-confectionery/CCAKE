# Usage of Singularity conservation_energy

1. Make sure that your system has singularity. In your personal computer, you
may simply install it. In a cluster, you need to load the respective module.
2. Using this folder as working directory, execute the `build_base_img.sh`
script to build the base image. This is an Ubuntu image with the required
packages installed.
3. Execute the `build.sh` script to build the image with the code. This will
clone the `kokkos` and `cabana` repositories and build them.
4. Execute the `run.sh` script to run the image. The code is not compiled yet,
but the image has all dependencies installed for an easy compilation.
5. Inside the container, navigate to `/CCake`, create a build directory and
execute `cmake ..` and `make -j`. This will compile the code.

> Note: The `/CCake` directory is mounted from the host machine. Any changes
> made to the code in the host machine is imediately reflected in the container.
