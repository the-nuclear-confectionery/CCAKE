#!/bin/bash

rm *.o
make
chmod +x charm.sh
cd df
rm *.o
make
chmod +x charm.sh
cd decays
rm *.o
make -f makefile reso
cd .. ; cd ..
