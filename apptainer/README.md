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
3. Execute the `build.sh` script to build the final image. It will create a new
   image with additional packages (at the moment, only TRENTo is provided) and
   execute the `bootstrap.sh` script in the root of the repository, which
   download the equation of state and compiles CCAKE.
4. Execute the `run.sh` script to run the image. You have access to a bash
   shell inside the container. The original repository in bound to the folder
   `/CCAKE` inside the container. A gubser test can be run as
   `mkdir -p output && ./ccake input/Input_Parameters_Gubser_checks.inp output`.

> Note: To build the image, one may need to have administrative privileges,
> depending on the configuration of the system. If for some reason you don't
> have these privileges, you can build the image in your local computer
> (steps 1 through 3) and copy the resulting CCAKE.sif
> file to the ./apptainer folder of the repository in the targeted cluster.
> Then you can proceed with steps 4 and 5.
>
> If you are going to perform a heavy computation (such as a full hydrodynamic
> simulation), it is recommended to run the image in a computing node.
> In UIUC CampusCluster, this command would be
> `srun -c 1 -n 1 -t 4:00:00 -p qgp --pty ./run.sh`. Similar commands can be
> used in other clusters that uses slurm as job scheduler.

## Development

The procedure above is optimized for using CCAKE. If you are modifying CCAKE,
you want to use the development base image, which contains extra tools (as gdb
for debugging and a ssh server to allow development tools like
`vscode` to use the container). Run the `build-base-img.sh`,
`build-base-img-dev.sh`, `build-dev.sh` and `run-dev.sh` in this sequence.

It is possible to use Visual Studio Code with these development containers. One
needs to add to their `~/.ssh/config` file
```
Host CCAKE~*
    RemoteCommand apptainer shell /storage/codes/BSQ/singularity/CCAKE-dev.sif
    RequestTTY yes

Host localhost CCAKE~localhost
    HostName 127.0.0.1
    User your-username
```

Install the remote development extensions on VSCode and configure it. To this end,
Add `remote.SSH.serverInstallPath` to the configuration the fields as below
```
"remote.SSH.serverInstallPath": {
  "CCAKE~localhost": "~/.vscode-container/CCAKE"
}
```
Also, make sure that `remote.SSH.enableRemoteCommand` and
`remote.SSH.useLocalServer` are enabled/set to `true`. You can now go
to under Remote Explorer in the Activity Bar and open `CCAKE-localhost`.
From there, you can open the cloned repository folder. If a password
is requested, it is the same as your user's password in the host system.

> Note: The `/CCAKE` directory is mounted from the host machine. Any changes
> made to the code in the host machine is immediately reflected in the container.
