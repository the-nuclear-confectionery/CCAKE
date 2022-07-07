#!/bin/bash

filename=$1

bash make_movie.sh $filename temperature
bash make_movie.sh $filename baryon_chemical_potential
bash make_movie.sh $filename strange_chemical_potential
bash make_movie.sh $filename electric_chemical_potential
bash make_movie.sh $filename energy_density
bash make_movie.sh $filename baryon_density
bash make_movie.sh $filename strange_density
bash make_movie.sh $filename electric_density
