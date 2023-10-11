#!/bin/bash

filename=$1

bash plotting_scripts/make_trajectory_plots.sh   $filename
bash plotting_scripts/make_density_plots.sh      $filename
bash plotting_scripts/make_eccentricity_plots.sh $filename
