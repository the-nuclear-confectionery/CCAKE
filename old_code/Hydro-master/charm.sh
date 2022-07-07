#!/bin/bash
<<COMMENT1
    run with:
    ./allnex.sh rnum 1stevent lastevent ideal/bulk/shear/shearbulk centrality 
    so,
    rnum=$1
    1stevent=$2
    lastevent=$3
    viscosity=$4
    centrality$5
    ic type (glb/cgc/gla/nex) $6
    rhic/lhc (if relevant) $7
COMMENT1
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
    v=""
    vv="i"
elif [ "$4" == "bulk" ] ; then
    v="bv"
    vv="bvc"
elif [  "$4" == "shearbulk" ] || [  "$4" == "shear" ]; then
    v="sbv"
    vv="sbvc"
fi


#echo First event $2 and final event $3
#echo outputfiles/"$t"/"$dec""$4"/"$5"/
for ((m="$2"; m<="$3"; m++))
do
    ./vusphydro input"$6""$1".dat "$1" "$m"
    sphpar=$(cat outputfiles/"$t"/"$dec""$4"/"$5"/temp_ic"$m".dat)
    sphout=$(wc -l < outputfiles/"$t"/"$dec""$4"/"$5"/"$v"freezeout_ev"$m".dat)
    if [ "$sphpar"-100 > "$sphout" ] ; then
       echo "$sphpar"  "$sphout"
       rm outputfiles/"$t"/"$dec""$4"/"$5"/"$v"freezeout_ev"$m".dat
       continue
    fi
    mv outputfiles/"$t"/"$dec""$4"/"$5"/"$v"freezeout_ev"$m".dat    df/input/"$t"/"$dec""$4"/"$5"/"$v"freezeout_ev"$m".dat
    cd df
    ./fo "$6"input"$1".dat "$1" "$m" 0 0 0
    ./decays/reso  out/"$t"/"$dec""$4"/"$5"/ev"$m""$vv"_dNdphidpp.dat out/"$t"/"$dec""$4"/"$5"/ev"$m"d"$vv"_dNdphidpp.dat reso16p.inp
    cd ..
done
