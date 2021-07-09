#!/bin/bash

maxFrames=250
resultsDirectory=outputfiles
physicalProfileStem=physical_quantities
suffix="_ev0.dat"
outsuffix="_ev0.png"

nFiles=$(\ls -1 ${resultsDirectory}/${physicalProfileStem}*suffix | wc -l)

for index in $(seq $nFiles)
do
	infile=${resultsDirectory}/${physicalProfileStem}${index}${suffix}

	# if file exists, generate corresponding frame(s)
	if [ -f "$infile" ] && [ $index -le $maxFrames ]
	then
		fwindex=`printf %03d $index`
		python3 plot_profile.py $infile ${resultsDirectory}/"eProfile"${index}${outsuffix} '$e$ (MeV/fm$^3$)' 7
		python3 plot_profile.py $infile ${resultsDirectory}/"rhoBProfile"${index}${outsuffix} '$\rho_B$ (fm$^{-3}$)' 8
		echo 'Plotted profiles from' $infile 'at' $(date)
	fi
done

framesPerSecond=25
pngs2mp4 $framesPerSecond ${resultsDirectory}/"eProfile"%03d${outsuffix} eProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoBProfile"%03d${outsuffix} rhoBProfile.mp4

