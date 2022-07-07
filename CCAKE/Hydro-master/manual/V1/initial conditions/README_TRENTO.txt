Generating initial conditions.

cd /projects/jn511_1/trento/build/src

1. Create a trento ".conf" file
2. mkdir #OUT_DIRECTOR_NAME#
3. Make a bash script #BASHSCRIPT#.sh that generates files
4. sbatch #BASHSCRIPT#.sh
5. From the folder /projects/jn511_1, run:
	mv trento/build/src/#OUT_DIRECTOR_NAME#/ v-USPhydro2/inputfiles/trento/
6. One now needs to generate the eccentricities (eccCM.dat) and initial entropy for each indidividual event (npart.dat) running eccs.sh.  
	./runECC.sh #OUT_DIRECTOR_NAME# #FILENAME# #END_FOLDER# #GRID_SPACING#
7. Combine all the different folders into one eccCM.dat and npart.dat.  In /projects/jn511_1
	./combine.sh trento/intputfiles/#OUT_DIRECTOR_NAME# eccCM.dat 0 #END_FOLDER#
	./combine.sh trento/intputfiles/#OUT_DIRECTOR_NAME# npart.dat 0 #END_FOLDER#
	cp trento/intputfiles/#OUT_DIRECTOR_NAME#/npart.dat v-USPhydro2/df/out/trento/#OUT_DIRECTOR_NAME#/


NOTES

1. Anything with #...# means you need to create your own name for it.
2. The cluster runs on slurm
3. To organize the initial conditions we separate them into subfolders titled 0, 1, 2 ...  with 999 initial conditions in each.  The entire code (and corresponding bash scripts) are set up to deal with this format. These initial conditions are NOT sorted in any special way.
4. Creating directors must be done in bash scripts, not within the code itself (the cluster is weird in that way). 
5. Trento bash scripts must include:
	module load mvapich2/2.1
	module load boost/1.60.0
for them to run properly.
6. The version of TRENTO on the cluster is NOT the same from the DUKE website.  Do NOT replace with the DUKE official version since the output has been updated to V-USPhydro format.  Additionally, we have upgraded the way some of the parameters for the different ions. 
7. One can change the different ion types on nucleus.cxx.  There one can play with beta2, beta3, beta4.  Make new ions with new names, don't delete old ones!
8.  It generally is very quick to generate all the initial conditions
9. npart.dat containts:
	event#   initial_entropy*
I know the name is misleading but it's what has been unfortunately hardcoded into the program a long time ago.  Also, initial_entropy is BEFORE multiplying by the constant, so the actual magnitude can be quite different. 
10. I assume that you're running the Equation of State with 2+1 quarks, so EOS21 is always used.  If you're playing around with a different type of EOS or viscosity or something that just call these folders something different. 
