#!/bin/bash

maxFrames=100
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
	${resultsDirectory}/"TProfile" \
	${insuffix} ${outsuffix} $nFrames '$T$ (MeV)' 3 temperature linear
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"muBProfile" \
	${insuffix} ${outsuffix} $nFrames '$\mu_B$ (MeV)' 4 baryon_chemical_potential linear
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"muSProfile" \
	${insuffix} ${outsuffix} $nFrames '$\mu_S$ (MeV)' 5 strange_chemical_potential linear
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"muQProfile" \
	${insuffix} ${outsuffix} $nFrames '$\mu_Q$ (MeV)' 6 electric_chemical_potential linear
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"eProfile" \
	${insuffix} ${outsuffix} $nFrames '$e$ (MeV/fm$^3$)' 7 energy_density log
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoBProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_B$ (fm$^{-3}$)' 8 baryon_density log
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoSProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_S$ (fm$^{-3}$)' 9 strange_density log
python3 plot_profile.py \
	${resultsDirectory}/${physicalProfileStem} \
	${resultsDirectory}/"rhoQProfile" \
	${insuffix} ${outsuffix} $nFrames '$\rho_Q$ (fm$^{-3}$)' 10 electric_density log


mkdir -p $moviesDirectory
rm -rf ${moviesDirectory}/eProfile.mp4 ${moviesDirectory}/rhoBProfile.mp4 \
       ${moviesDirectory}/rhoSProfile.mp4  ${moviesDirectory}/rhoQProfile.mp4 \
       ${moviesDirectory}/TProfile.mp4 ${moviesDirectory}/muBProfile.mp4 \
       ${moviesDirectory}/muSProfile.mp4  ${moviesDirectory}/muQProfile.mp4

rm -f ${moviesDirectory}/TProfile.mp4
       
pngs2mp4 $framesPerSecond ${resultsDirectory}/"TProfile"%03d${outsuffix} ${moviesDirectory}/TProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"muBProfile"%03d${outsuffix} ${moviesDirectory}/muBProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"muSProfile"%03d${outsuffix} ${moviesDirectory}/muSProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"muQProfile"%03d${outsuffix} ${moviesDirectory}/muQProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"eProfile"%03d${outsuffix} ${moviesDirectory}/eProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoBProfile"%03d${outsuffix} ${moviesDirectory}/rhoBProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoSProfile"%03d${outsuffix} ${moviesDirectory}/rhoSProfile.mp4
pngs2mp4 $framesPerSecond ${resultsDirectory}/"rhoQProfile"%03d${outsuffix} ${moviesDirectory}/rhoQProfile.mp4

