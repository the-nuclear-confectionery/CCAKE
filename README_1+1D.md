Please make sure to sync with the most recent version of the 'skv' branch because I have made the following changes:
- RK 4 has been implemented and can be run by setting rk_order in the yaml file to 4, if desired.
- Energy conservation has been implemented for rapidity and each time step will show the conservation status as the simulation runs.
- Added a missing tau in grad_kernel calculation in kernel.cpp.

To run 1+1 D:
- To create 1+1 D simulation, the initial conditions needs to be created using long_ic_generator.py file in the CCAKE/utilities folder using the command $ python3 long_ic_generator.py 1.0, where 1.0 is the tau for that solution. That value can be changed to create solutions for other taus. This initial condition file will be saved in the folder /CCAKE/long_data. If this folder doesn't presently exist, it will need to be created first. For Bjorken, the steps are the same, except the initial_conditions will be saved in the folder /CCAKE/bjorken and the IC generator file is name bjorken_ic_generator.py.

- The configuration files for both soultions are in the /CCAKE/config folder with the names config_long.yaml (for 1+1 D) and config_bjorken.yaml. Make sure to set the runge-kutta order you would like to use and change hT if you want to change the smoothing parameter. heta is not being used yet and hT is being treated as heta for 1+1 D purposes. Will be fixed when it is time to run 3+1 D.
  
- To run the 1+1 D simulation, use the command $ srun -A bbkr-delta-cpu --time=01:00:00 --cpus-per-task=4 --partition=cpu --pty make run-long. The results will be saved in the folder /CCAKE/build/output_long.

- To run Bjorken, use the command srun -A bbkr-delta-cpu --time=01:00:00 --cpus-per-task=4 --partition=cpu --pty make run-bjorken. The results will be saved in the folder /CCAKE/build/output_bjorken.

- The plotting scripts to generate energy density vs eta and u_eta vs eta plots have been added to the /CCAKE/plotting_scripts folder.
  - plot_long.py is suited best to run on the cluster and will generate energy density vs eta and u_eta vs eta plots       for different taus. To run it, the command $ plot_long.py <path to analytical solution folder> <path to sim 
    solution folder> can be used.
  The following three scrips are direct uploads from jupyter notebook, but can be referred to in the case that           plot_long.py needs to be adapted to generate comparison plots for, say, smoothing parameter or runge-kutta order. It   would be worth going through simple plots.py to improve the formatting and labelling in plot_long.py, but it is not    entirely necessary.
  - simple plots.py will generate energy density vs eta and u_eta vs eta for different taus.
  - RK Comparison.py will generate plots over different taus to compare results for the other runge-kutta orders being     used in the code.
  - comparative plots.py generates plots for different smoothing values.


