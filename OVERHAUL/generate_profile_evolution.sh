#!/bin/bash

maxFrames=300
framesPerSecond=30

resultsDirectory=results
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
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoSProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_S$ (fm$^{-3}$)' 9
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoQProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_Q$ (fm$^{-3}$)' 10

mkdir -p $moviesDirectory
rm -rf ${moviesDirectory}/eProfile.mp4 ${moviesDirectory}/rhoBProfile.mp4 \
       ${moviesDirectory}/rhoSProfile.mp4  ${moviesDirectory}/rhoQProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"eProfile"%03d${outsuffix} ${moviesDirectory}/eProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoBProfile"%03d${outsuffix} ${moviesDirectory}/rhoBProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoSProfile"%03d${outsuffix} ${moviesDirectory}/rhoSProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoQProfile"%03d${outsuffix} ${moviesDirectory}/rhoQProfile.mp4

