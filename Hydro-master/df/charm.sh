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

#make
#cd decays 
#make -f makefile reso 
#cd .. 

if [ "$6" == "glb" ] ; then
	t="glauber"
elif [ "$6" == "cgc" ] ; then
	t="cgc"
elif [ "$6" == "nex" ] ; then
	t="nexus"
elif [ "$6" == "gla" ] ; then
	t="glasma"
elif [ "$6" == "trento" ] ; then
        t="trento"
fi
if [ "$7" != "" ] ; then
	dec="$7"/
else 
	dec=""
fi
if [ "$4" == "ideal" ] ; then
	v="i"
elif [ "$4" == "bulk" ] ; then
	v="bvc"
elif [  "$4" == "shearbulk" ] || [  "$4" == "shear" ]; then
	v="sbvc"
fi
echo First event $2 and final event $3
for  (( k=$2; k<=$3; k++ ))
do
    echo "$k"
./fo "$6"input"$1".dat "$1" "$k" 0 0 0
./decays/reso  out/"$t"/"$dec""$4"/"$5"/ev"$k""$v"_dNdphidpp.dat out/"$t"/"$dec""$4"/"$5"/ev"$k"d"$v"_dNdphidpp.dat reso16p.inp
done
