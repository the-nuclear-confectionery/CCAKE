#!/bin/bash

maxFrames=250
framesPerSecond=25

resultsDirectory=outputfiles
physicalProfileStem=physical_quantities
insuffix="_ev0.dat"
outsuffix="_ev0.png"
moviesDirectory=movies

nFiles=$(\ls -1 ${resultsDirectory}/${physicalProfileStem}*${insuffix} | wc -l)

# old version which looped over each file one at a time
#for index in $(seq $nFiles)
#do
#	infile=${resultsDirectory}/${physicalProfileStem}${index}${insuffix}
#
#	# if file exists, generate corresponding frame(s)
#	if [ -f "$infile" ] && [ $index -le $maxFrames ]
#	then
#		fwindex=`printf %03d $index`
#		python3 plot_profile.py $infile ${resultsDirectory}/"eProfile"${index}${outsuffix} '$e$ (MeV/fm$^3$)' 7
#		python3 plot_profile.py $infile ${resultsDirectory}/"rhoBProfile"${index}${outsuffix} '$\rho_B$ (fm$^{-3}$)' 8
#		echo 'Plotted profiles from' $infile 'at' $(date)
#	fi
#done

# use the minimum of these two quantities
nFrames=$((nFiles<maxFrames ? nFiles : maxFrames))

python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"eProfile" \
	${insuffix} ${outsuffix} $nFrames '$e$ (MeV/fm$^3$)' 7
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoBProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_B$ (fm$^{-3}$)' 8



mkdir -p $moviesDirectory
pngs2mp4 $framesPerSecond ${resultsDirectory}/"eProfile"%03d${outsuffix} ${moviesDirectory}/eProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoBProfile"%03d${outsuffix} ${moviesDirectory}/rhoBProfile.mp4

