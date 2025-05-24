#!/bin/bash
#----------

BUILD_DIR=build

mkdir -p ${BUILD_DIR}
( cd ${BUILD_DIR} \
  && cmake -DCMAKE_BUILD_TYPE=Debug ..     \
  && make         \
  && make install )
