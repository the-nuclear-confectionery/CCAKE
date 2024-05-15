# CCAKE

CCAKE is a relativistic viscous hydrodynamic code with 3 conserved charges 
(baryon number, strangeness, and electric charge) that uses Smoothed Particle 
Hydrodynamics. CCAKE can make state-of-the-art predictions for heavy-ion 
collisions. It uses the 4D lattice QCD equation of state with T, baryon, 
strangeness, and electric charge that is coupled to the PDG16+ particle
list.

If you use this code in your research, please remember to cite us.

## 1. Apptainer (Singularity) usage

It is recommended to use the apptainer (former singularity) to execute the code. 
You can find out instructions in [here](apptainer/README.md).

## 2. Compiling from source instructions

### 2.1. Dependencies

To compile `ccake`, you need to have the following dependencies installed:

- C++ compiler with support for C++17 (g++ > 7.2.0)
- CMake (>= 3.4)
- HDF5
- GSL

> Dealing with missing dependencies: The easiest way is to install the missing
> dependencies using your package manager. For example, on Ubuntu, you can run
> `sudo apt install cmake libhdf5-dev libgsl-dev cmake`. If for some reason you
> you can`t or don't want to install the dependencies (e.g., you are on a
> cluster and don't have root access), you can use either compile the
> dependencies from source or use a virtual environment (recommended). HPC
> clusters may also offer them as loadable modules.
>
> To use the virtual environment, you need to have
> [`mamba`](https://mamba.readthedocs.io/en/latest/) installed. Then, run
> `mamba env create -f environment.yml` to create the environment. To activate
> it, run `conda activate ccake`. To deactivate it, run `conda deactivate`.

### 2.2. Compiling

The recommended way to compile is by using the `build_all.sh` script. You can
manually compile the project by following the steps below.

1. Create a build directory and enter it: `mkdir build && cd build`
2. Configure the build with `cmake ..`
3. Build with `make -jN` where `N` is the number of cores you want to use. For 
   a personal computer, `N=4` is a good choice.

## 3. Usage

Make sure that in your working directory, a copy or a symlink of the folder
`EoS` is present. Create an output directory and then run `ccake` as
`ccake path/to/input-file path/to/output-directory`. For example, if you are
in the project root directory, you can run
`build/bin/ccake input/input/Input_Parameters_ccake.inp output`.

## 4. Reporting bugs and getting help

If you are facing some problem or found some bug, feel free to open an 
[issue](https://github.com/the-nuclear-confectionery/CCAKE/issues) or create
a pull request.

## 5. Debugging

If you are developing `ccake`, you might want to debug it. To do so, you need
to compile it with debug symbols. Steps to do so are as follows:

1. Create a build directory and enter it: `mkdir build && cd build` as above
2. At the configure step, add the flag `-DCMAKE_BUILD_TYPE=Debug` to the
   `cmake` command. For example, `cmake -DCMAKE_BUILD_TYPE=Debug ..`
3. Build with `make -jN` as above

To debug, you can use `gdb` or `lldb`. For example,
run `gdb --args build/bin/ccake path/to/input-file`. If it is your first time
debugging, you might want to read the
[gdb tutorial](https://sourceware.org/gdb/onlinedocs/gdb/).

> VSCode users: You can use the `C/C++` extension to debug `ccake`. Install it
> from the extensions tab. A launch configuration is already provided in the
> `.vscode` directory. You can set breakpoints by clicking on the left of the
> line number (a red dot will appear). Then, run the debugger by pressing `F5`.
> You may inspect the variables by hovering over them or in the left panel.

## 6. Generating doxygen documentation

To generate the doxygen documentation, you need to have doxygen installed. Then,
in the project root directory, run `doxygen Doxyfile`. The documentation will
be generated in the `docs` directory. Two formats for the documentation will be
generated: html and latex. To view the html documentation, open the file
`docs/html/index.html` in your browser. To view the latex documentation, you
need to compile the latex files. To do so, run `make` in the `docs/latex`. The
generated pdf will be in `docs/latex/refman.pdf`.
