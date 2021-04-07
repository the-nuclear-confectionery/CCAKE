#!/bin/bash 

./vns/fo "$1" "$2" "$3" "$4"  "$5" "$6"  "$7"
cd out
c++ main.cpp
cd ..
if [ "$7" == 3 ] ; then
	let t=0
else
	let t=1
fi



./out/a.out "$1" "$2" "$3" "$4"   "$t"  "$7"


