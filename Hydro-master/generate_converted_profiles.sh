#!/bin/bash

resultsDirectory=outputfiles
SPHprofileStem=bsqsveprofile
physicalProfileStem=physical_quantities
suffix="_ev0.dat"

nFiles=$(\ls -1 ${resultsDirectory}/${SPHprofileStem}*${suffix} | wc -l)

# compile from scratch since it's so quick
rm ./convert_SPH_to_physical
g++ -std=c++11 convert_SPH_to_physical.cpp -o convert_SPH_to_physical

for index in $(seq $nFiles)
do
	if [ ! -f "${resultsDirectory}/${physicalProfileStem}${index}${suffix}" ]
	then
		./convert_SPH_to_physical \
			$SPHprofileStem $physicalProfileStem $index
		echo 'Converted' $index 'of' $nFiles 'profiles at' $(date)
	fi
done
