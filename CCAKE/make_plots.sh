#!/bin/bash

filename=$1

bash make_trajectory_plots.sh   $filename
bash make_density_plots.sh      $filename
bash make_eccentricity_plots.sh $filename
