#!/bin/bash
<<COMMENT1
    run with:
    ./ebedf.sh rnum 1stevent lastevent ideal/bulk/shear/shearbulk centrality 
    so,
    rnum=$1
    1stevent=$2
    lastevent=$3
    viscosity=$4
    centrality$5
    ic type $6
    rhic/lhc (if relevant) $7
COMMENT1
if [ "$4" == "ideal" ] ; then
	v="i"
elif [ "$4" == "bulk" ] ; then
	v="bvc"
elif [  "$4" == "shearbulk" ] || [  "$4" == "shear" ]; then
	v="sbvc"
fi

if [ "$6" == "glb" ] ; then
	t="glauber"
elif [ "$6" == "cgc" ] ; then
	t="cgc"
elif [ "$6" == "nex" ] ; then
	t="nexus"
elif [ "$6" == "gla" ] ; then
	t="glasma"
fi

if [ "$7" != "" ] ; then
	dec="$7"/
else 
	dec=""
fi

cd decays 
make -f makefile reso 
cd .. 
for (( j=$2; j<=$3; j++ ))
do
   ./decays/reso  out/"$t"/"$dec""$4"/"$5"/ev"$j""$v"_dNdphidpp.dat out/"$t"/"$dec""$4"/"$5"/ev"$j"d"$v"_dNdphidpp.dat 1 1
done 
