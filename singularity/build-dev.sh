#!/usr/bin/env bash
apptainer build --fakeroot --sandbox ccake-dev.sif ccake-dev.def
#singularity instance start ccake.sif ccake