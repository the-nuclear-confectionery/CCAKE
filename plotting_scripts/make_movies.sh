#!/bin/bash

filename=$1

bash plotting_scripts/make_movie.sh $filename temperature
bash plotting_scripts/make_movie.sh $filename baryon_chemical_potential
bash plotting_scripts/make_movie.sh $filename strange_chemical_potential
bash plotting_scripts/make_movie.sh $filename electric_chemical_potential
bash plotting_scripts/make_movie.sh $filename energy_density
bash plotting_scripts/make_movie.sh $filename baryon_density
bash plotting_scripts/make_movie.sh $filename strange_density
bash plotting_scripts/make_movie.sh $filename electric_density
