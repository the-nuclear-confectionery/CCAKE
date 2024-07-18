#!/bin/bash
#----------

BUILD_DIR=build

mkdir -p ${BUILD_DIR}
( cd ${BUILD_DIR} \
  && cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang ..     \
  && make         \
  && make install )