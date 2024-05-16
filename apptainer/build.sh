#!/usr/bin/env bash
singularity build --fakeroot ccake.sif ccake.def
singularity run --fakeroot --cleanenv --no-home -B ../:/CCAKE ccake.sif \
bash -c "cd /CCAKE && source bootstrap.sh"