#!/usr/bin/env bash
singularity build --fakeroot ubuntu.sif ubuntu.def
#singularity instance start ccake.sif ccake