#!/usr/bin/env bash
singularity build --fakeroot ccake.sif ccake.def
#singularity instance start ccake.sif ccake