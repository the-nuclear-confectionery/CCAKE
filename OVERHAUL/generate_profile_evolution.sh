#!/bin/bash

maxFrames=250
framesPerSecond=2

resultsDirectory=outputfiles
physicalProfileStem=physical_quantities
insuffix=".dat"
outsuffix=".png"
moviesDirectory=movies

nFiles=$(\ls -1 ${resultsDirectory}/${physicalProfileStem}*${insuffix} | wc -l)

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
rm -rf ${moviesDirectory}/eProfile.mp4 ${moviesDirectory}/rhoBProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"eProfile"%03d${outsuffix} ${moviesDirectory}/eProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoBProfile"%03d${outsuffix} ${moviesDirectory}/rhoBProfile.mp4

