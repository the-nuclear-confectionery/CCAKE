#!/usr/bin/env bash
singularity build --fakeroot ccake.sif ../defs/ccake.def
#singularity instance start ccake.sif ccake
