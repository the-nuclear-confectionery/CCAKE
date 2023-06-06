# Usage of Apptainer

> Note: Apptainer is just a new name for singularity. The name change was made 
> after they joined the Linux Foundation at the end of November, 2021. Many 
> systems still uses the singularity name. You can safely replace `apptainer`
> by `singularity` in these systems. An easy way of doing so is by creating 
> an alias in your bashrc file.

1. Make sure that your system has `apptainer` (or singularity). In your personal 
   computer, you  may simply install it. In a cluster, you need to load the 
   respective module.
2. Using this folder as working directory, execute the `build-base-img.sh`
   script to build the base image. This is an Ubuntu image with the required
   packages installed.
3. Execute the `build.sh` script to build the image with the code. This will
   clone the `kokkos` and `cabana` repositories and build them.
4. Execute the `run.sh` script to run the image. The code is not compiled yet,
   but the image has all dependencies installed for an easy compilation.
5. Inside the container, navigate to `/CCake`, create a build directory and
   execute `cmake ..` and `make -j`. This will compile the code.

> Note: The `/CCake` directory is mounted from the host machine. Any changes
> made to the code in the host machine is immediately reflected in the container.

## Development

The procedure above is optimized for using CCake. If you are modifying CCake,
you want to use the development base image, which contains extra tools (as 
a ssh server to allow development tools like `vscode` to use the container and
gdb). Run the `build-base-img.sh`, `build-base-img-dev.sh`, `build-dev.sh` and
`run-dev.sh` in this sequence.

It is possible to use Visual Studio Code with these development containers. One 
needs to add to their `~/.ssh/config` file
```
Host ccake~*
    RemoteCommand apptainer shell /storage/codes/BSQ/singularity/ccake-dev.sif
    RequestTTY yes

Host localhost ccake~localhost
    HostName 127.0.0.1
    User your-username
```

Install the remote development extensions on VSCode and configure it. To this end,
Adding to `remote.SSH.serverInstallPath`  configuration the fields as below
```
"remote.SSH.serverInstallPath": {
  "ccake~localhost": "~/.vscode-container/ccake"
}
```
Also, make sure that `remote.SSH.enableRemoteCommand` and 
`remote.SSH.useLocalServer` are enabled/set to `true`. You can now go
to under Remote Explorer in the Activity Bar and open `ccake-localhost`.
From there, you can open the cloned repository folder. If a password
is requested, it is the same as your user's password in the host system.