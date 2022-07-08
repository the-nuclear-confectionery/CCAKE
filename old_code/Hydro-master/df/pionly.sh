#!/bin/bash 

for  (( k=$1; k<=$2; k++ ))
do
./fo trentoinputpPb5.dat pPb5 "$k" 0 0
done
