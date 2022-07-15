#!/bin/bash
#----------

BUILD_DIR=build

mkdir -p ${BUILD_DIR}
( cd ${BUILD_DIR} \
  && cmake ..     \
  && make VERBOSE=1;        \
  && make install )